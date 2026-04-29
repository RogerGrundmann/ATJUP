/*
 * Atmosphere General Circulation Model (ATJUP)
 * Finite difference scheme for the 3D Navier-Stokes equations on Jupiter
 * 4th-order Runge-Kutta integration
 *
 * Initial and boundary conditions for thermodynamic and chemical species fields
 */
#include "cJupiterModel.h"
#include "Utils.h"
#include "SaturationAdjustmentJup.h"

using namespace std;
using namespace JupiterUtils;

// ============================================================================
// Physical and Numerical Constants
// ============================================================================
namespace JupiterInitConstants {
    // H2S density correction to match Planetary Sciences p. 90 (2010)
    constexpr double CORR_H2S   = 4.0e-3;
    // NH4SH density correction to match Planetary Sciences p. 90 (2010)
    constexpr double CORR_NH4SH = 6.0;
    // Saturation vapour scaling factor (empirical)
    constexpr double MAGNUS_COEFF = 2.0;
    // Pole-to-equator r_max ratio
    constexpr double R_MAX_POLE_FRAC = 0.7;
    // NH4SH onset temperature tolerance
    constexpr double NH4SH_TEMP_MARGIN = 1.05;
}


using namespace JupiterInitConstants;
using namespace JupiterUtils;


    cJupiterModel SaturationAdjustmentJup;



void cJupiterModel::init_tropopause_layers(){                                                                                                                                                         
    cout << endl << endl << endl << "      Jupiter: init_tropopause_layers" << endl;                                                                                                                        
                                                                                                                                                                                                           

    int j_max = jm - 1;                                                                                                                                                                                  
    int j_half = j_max / 2;                                                                                                                                                                              
                  
    // Derive x_max so that Agnesi(tropopause_equator, x_max) == tropopause_pole exactly.                                                                                                                
    // Agnesi: a^3/(a^2+x^2) = b  =>  x = a * sqrt(a/b - 1)
    // Requires tropopause_equator > tropopause_pole (always true physically).                                                                                                                           
    double x_max = tropopause_equator                                                                                                                                                                    
                   * std::sqrt(tropopause_equator / tropopause_pole - 1.0);                                                                                                                              
                  
  cout << "tropopause_pole=" << tropopause_pole                                                                                                                                                            
       << " x_max=" << x_max                                                                                                                                                                               
       << " pole_index=" << round(tropopause_pole/L_atm) << endl;                                                                                                                                          



                                                                                                                                                                                         
    // Build symmetric cache of heights [m] and grid indices in one pass.                                                                                                                                
    std::vector<double> tropo_height_cache(jm);
    tropopause_layers = std::vector<double>(jm);                                                                                                                                                         
  
    for(int j = 0; j <= j_half; j++){                                                                                                                                                                    
        double x = x_max * (double)(j_half - j) / (double)j_half;
        double h = JupiterUtils::Agnesi(tropopause_equator, x);                                                                                                                                             
        tropo_height_cache[j]       = h;
        tropo_height_cache[j_max-j] = h;                                                                                                                                                               
        tropopause_layers[j]        = round(h / L_atm);                                                                                                                                                 
        tropopause_layers[j_max-j]  = tropopause_layers[j];                                                                                                                                             
    }                                                                                                                                                                                                    
/*                                                                                                                                                       
    #pragma omp parallel for schedule(static)                                                                                                                                                            
   for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){                                                                                                                                                                     
            Tropopause.y[j][k] = tropo_height_cache[j];
        }                                                                                                                                                                                                
    }           
*/                                                                                                                                                                                             
    cout << "      Jupiter: init_tropopause_layers ended" << endl;                                                                                                                                          
}





// ============================================================================
// Temperature Initialization
// ============================================================================
void cJupiterModel::init_temperature(){
    std::cout << "\n\n\n      ATJUP: init_temperature" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const int    j_half    = (jm - 1) / 2;
    const double d_j_half  = (double)j_half;
    const double t_h2_eff  = t_pole - t_equator;

    // ========================================================================
    // Latitudinal surface temperature + vertical lapse-rate profile
    // ========================================================================
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
//            const double d_j = (double)j;
//            double ratio = (double)j / d_j_half;
            // Parabolic pole-to-equator distribution (non-dimensional)
//            t.x[0][j][k] = (t_h2_eff * (d_j * d_j / (d_j_half * d_j_half)
//                            - 2.0 * d_j / d_j_half) + t_pole) / t_ref;
            double ratio = (double)j / d_j_half;
            t.x[0][j][k] = (t_h2_eff * JupiterUtils::parabola(ratio) + t_pole) / t_ref;

            // Linear lapse-rate profile upward from surface
            for (int i = 1; i < im; i++) {
                const double height = get_layer_height(i);
                t.x[i][j][k] = -gam * height / t_ref + t.x[0][j][k];
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for init_temperature\n", elapsed.count() * 1e-9);
    std::cout << "      ATJUP: init_temperature ended" << std::endl;
}


// ============================================================================
// Dynamic Pressure Initialization (perturbation pressure — starts at zero)
// ============================================================================
void cJupiterModel::init_PressureDynamic(){
    std::cout << "\n\n\n      ATJUP: init_PressureDynamic" << std::endl;

    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            for (int i = 0; i < im; i++) {
                p_dyn.x[i][j][k] = 0.0;
            }
        }
    }

    std::cout << "      ATJUP: init_PressureDynamic ended" << std::endl;
}


// ============================================================================
// Static Pressure Initialization (hypsometric / power-law profile)
// ============================================================================
void cJupiterModel::init_PressureStatic(){
    std::cout << "\n\n\n      ATJUP: init_PressureStatic" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double exp_pressure = g / (gam * R_ref);

    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            for (int i = 0; i < im; i++) {
                const double t_u = t.x[i][j][k] * t_ref;
                p_stat.x[i][j][k] = p_ref * pow(t_u / t_ref, exp_pressure);
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for init_PressureStatic\n", elapsed.count() * 1e-9);
    std::cout << "      ATJUP: init_PressureStatic ended" << std::endl;
}


// ============================================================================
// Body-Force Fields (Coriolis, Centrifugal, Buoyancy, Pressure-Gradient)
// ============================================================================
void cJupiterModel::Forces(){
    std::cout << "\n\n\n      ATJUP: Forces" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Interior points only (boundaries handled elsewhere)
    // ========================================================================
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 1; i < im-1; i++) {
        for (int j = 1; j < jm-1; j++) {
            const double rm       = rad.z[i];
            double sinthe         = sin(the.z[j]);
            double costhe         = cos(the.z[j]);
            if (j > 90) costhe    = -costhe;
            const double rmsinthe = rm * sinthe;

            for (int k = 1; k < km-1; k++) {

                // Coriolis components
                const double Cor_r   = -2.0 * omega * sinthe * w.x[i][j][k];
                const double Cor_the = +2.0 * omega * costhe * w.x[i][j][k];
                const double Cor_phi = +2.0 * omega * (-costhe * v.x[i][j][k]
                                        + sinthe * u.x[i][j][k]);

                // Pressure-gradient components (central differences)
                const double dpdr   = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) / (2.0 * dr);
                const double dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) / (2.0 * dthe);
                const double dpdphi = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) / (2.0 * dphi);

                CoriolisForce.x[i][j][k] = Coriolis * r_mix
                    * sqrt((pow(Cor_r, 2) + pow(Cor_the, 2) + pow(Cor_phi, 2)) / 3.0);

                CentrifugalForce.x[i][j][k] = centrifugal * r_mix
                    * omega * omega * rm * (1.0 + fabs(sinthe));

                BuoyancyForce.x[i][j][k] = buoyancy
                    * r_mix * g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])
                    / (r_mix * R_mix * t.x[i][j][k] * t_ref) * 1e5;

                PresGradForce.x[i][j][k] =
                    -sqrt((pow(dpdr, 2)
                          + pow(dpdthe / rm, 2)
                          + pow(dpdphi / rmsinthe, 2)) / 3.0) / L_atm * 1.0e5;
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for Forces\n", elapsed.count() * 1e-9);
    std::cout << "      ATJUP: Forces ended" << std::endl;
}


// ============================================================================
// Vapour, Cloud and Ice Initialization (H2O and NH3)
// ============================================================================
void cJupiterModel::init_vapour_cloud_ice(std::string gas,
    double &c_tropopause, double &coeff_A, double &coeff_B,
    double &coeff_A_i, double &coeff_B_i,
    double &t_0, double &t_00,
    double &ep, double &r, double &m,
    double &C, double &L0, double &R,
    double &del_alf, double &del_bet,
    Array &c, Array &cloud, Array &ice){


    std::cout << "\n\n\n      ATJUP: init_vapour_cloud_ice of " << gas << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Per-latitude distribution parameters
    // ========================================================================
    double r_max_equator       = 0.0;
    double r_max_pole          = 0.0;
    double r_max_add_equator   = 0.0;
    double r_max_add_pole      = 0.0;
    double t_add_equator       = 0.0;
    double t_add_pole          = 0.0;
    double cloud_loc_equator   = 0.0;
    double cloud_loc_pole      = 0.0;
    double magnus              = MAGNUS_COEFF;

    if (gas == "H2O") {
        r_max_equator      = r_h2o;
        r_max_pole         = R_MAX_POLE_FRAC * r_h2o;
        cloud_loc_equator  = 1.0;
        cloud_loc_pole     = 0.0;
    }
    if (gas == "NH3") {
        r_max_equator      = r_nh3;
        r_max_pole         = R_MAX_POLE_FRAC * r_nh3;
        r_max_add_equator  = r_nh3_add;
        r_max_add_pole     = R_MAX_POLE_FRAC * r_nh3_add;
        t_add_equator      = t_0_nh3;
        t_add_pole         = 1.01 * t_add_equator;
        cloud_loc_equator  = 9.0;
        cloud_loc_pole     = 7.0;
    }

    // Build latitude-dependent lookup vectors (sequential: jm is small)
    const double d_j_half      = (double)(jm - 1) / 2.0;
    const double cloud_loc_eff = cloud_loc_pole - cloud_loc_equator;
    const double r_max_eff     = r_max_pole - r_max_equator;
    const double r_max_add_eff = r_max_add_pole - r_max_add_equator;
    const double t_add_eff     = t_add_pole - t_add_equator;

    cloud_loc  = std::vector<double>(jm, cloud_loc_pole);
    r_max      = std::vector<double>(jm, r_max_pole);
    r_max_add  = std::vector<double>(jm, r_max_add_pole);
    t_add      = std::vector<double>(jm, t_add_pole);

    for (int j = 0; j < jm; j++) {
        const double ratio = (double)j / d_j_half;
        const double par   = JupiterUtils::parabola(ratio);
        cloud_loc[j]  = cloud_loc_eff  * par + cloud_loc_pole;
        r_max[j]      = r_max_eff      * par + r_max_pole;
        r_max_add[j]  = r_max_add_eff  * par + r_max_add_pole;
        t_add[j]      = t_add_eff      * par + t_add_pole;
    }

    // ========================================================================
    // Main initialisation: vapour field
    // ========================================================================
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {

            for (int i = 0; i <= im-1; i++) {
                double t_u          = t.x[i][j][k] * t_ref;
                const double p_u    = p_stat.x[i][j][k];
                const double E_Rain = SaturationAdjustmentJup::saturation_vapour_pressure(
                                          t_u, C, L0, R, del_alf, del_bet);
                const double q_Rain = ep * E_Rain / p_u;

                double cv = magnus * r_mix * q_Rain;
                if (cv >= r_max[j])   cv = r_max[j];
                if ((gas == "NH3") && (t_u > t_add[j]))  cv = r_max_add[j];
                c.x[i][j][k] = cv;
            }
        }
    }

    // ========================================================================
    // NH3 smoothing at warm-layer transition (Gauss-Seidel stencil — serial)
    // ========================================================================
    for (int j = 2; j < jm-2; j++) {
        for (int k = 2; k < km-2; k++) {
            for (int i = 2; i < im-2; i++) {
                const double t_u = t.x[i][j][k] * t_ref;
                if ((gas == "NH3") && (t_u >= t_add[j])) {
                    c.x[i][j][k] =
                        (c.x[i+1][j][k] + c.x[i-1][j][k]
                       + c.x[i+2][j][k] + c.x[i-2][j][k]
                       + c.x[i][j+1][k] + c.x[i][j-1][k]
                       + c.x[i][j+2][k] + c.x[i][j-2][k]
                       + c.x[i][j][k+1] + c.x[i][j][k-1]
                       + c.x[i][j][k+2] + c.x[i][j][k-2]) / 12.0;
                }
            }
        }
    }

    // ========================================================================
    // Zero all species inside solid terrain
    // ========================================================================
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            for (int i = 0; i < im; i++) {
                if (is_land(SeaMount, i, j, k)) {
                    c.x[i][j][k]          = 0.0;
                    cloud.x[i][j][k]      = 0.0;
                    ice.x[i][j][k]        = 0.0;
                }
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for init_vapour_cloud_ice (%s)\n",
           elapsed.count() * 1e-9, gas.c_str());
    std::cout << "      ATJUP: init_vapour_cloud_ice of " << gas << " ended" << std::endl;
}


// ============================================================================
// H2S Vapour Initialization
// ============================================================================
void cJupiterModel::init_h2s(std::string gas,
    double &c_tropopause, double &coeff_A, double &coeff_B,
    double &coeff_A_i, double &coeff_B_i,
    double &t_0, double &t_00,
    double &ep, double &r, double &m,
    double &C, double &L0, double &R,
    double &del_alf, double &del_bet, Array &c){


    std::cout << "\n\n\n      ATJUP: init_h2s of " << gas << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Latitude-dependent density cap
    // ========================================================================
    const double r_max_equator = r_h2s;
    const double r_max_pole    = R_MAX_POLE_FRAC * r_h2s;
    const double d_j_half      = (double)(jm - 1) / 2.0;
    const double r_max_eff     = r_max_pole - r_max_equator;

    r_max = std::vector<double>(jm, r_max_pole);
    for (int j = 0; j < jm; j++) {
        r_max[j] = r_max_eff * JupiterUtils::parabola((double)j / d_j_half) + r_max_pole;
    }

    // ========================================================================
    // Vapour field — independent columns, safe to parallelise
    // ========================================================================
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            for (int i = 0; i <= im-1; i++) {
                double t_u          = t.x[i][j][k] * t_ref;
                const double p_u    = p_stat.x[i][j][k];
                const double E_Rain = SaturationAdjustmentJup::saturation_vapour_pressure(
                                          t_u, C, L0, R, del_alf, del_bet);
                const double q_Rain = ep * E_Rain / p_u;

                double cv = CORR_H2S * r_mix * q_Rain;
                if (cv >= r_max[j])  cv = r_max[j];
                c.x[i][j][k] = cv;
            }
        }
    }

    // ========================================================================
    // Zero inside terrain
    // ========================================================================
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            for (int i = 0; i < im; i++) {
                if (is_land(SeaMount, i, j, k))  c.x[i][j][k] = 0.0;
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for init_h2s (%s)\n",
           elapsed.count() * 1e-9, gas.c_str());
    std::cout << "      ATJUP: init_h2s of " << gas << " ended" << std::endl;
}


// ============================================================================
// NH4SH Vapour Initialization
// ============================================================================
void cJupiterModel::init_nh4sh(std::string gas, double &c_tropopause,
    double &coeff_A, double &coeff_B,
    double &ep, double &r, double &m,
    double &C, double &L0, double &R,
    double &del_alf, double &del_bet, Array &c){


    std::cout << "\n\n\n      ATJUP: init_nh4sh of " << gas << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Latitude-dependent density cap
    // ========================================================================
    const double r_max_equator = r_nh4sh;
    const double r_max_pole    = R_MAX_POLE_FRAC * r_nh4sh;
    const double d_j_half      = (double)(jm - 1) / 2.0;
    const double r_max_eff     = r_max_pole - r_max_equator;

    r_max = std::vector<double>(jm, r_max_pole);
    for (int j = 0; j < jm; j++) {
        r_max[j] = r_max_eff * JupiterUtils::parabola((double)j / d_j_half) + r_max_pole;
    }

    // ========================================================================
    // Vapour field — independent columns, safe to parallelise
    // ========================================================================
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            for (int i = 0; i < im; i++) {
                double t_u          = t.x[i][j][k] * t_ref;
                const double p_u    = p_stat.x[i][j][k];
                const double E_Rain = SaturationAdjustmentJup::saturation_vapour_pressure(
                                          t_u, C, L0, R, del_alf, del_bet);
                const double q_Rain = ep * E_Rain / p_u;

                double cv = CORR_NH4SH * r_mix * q_Rain;
                if (cv <= 0.0)                       cv = 0.0;
                if (cv >= r_max[j])                  cv = r_max[j];
                if (t_u > NH4SH_TEMP_MARGIN * t_0_nh4sh)  cv = 0.0;
                if (is_land(SeaMount, i, j, k))      cv = 0.0;
                c.x[i][j][k] = cv;
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for init_nh4sh (%s)\n",
           elapsed.count() * 1e-9, gas.c_str());
    std::cout << "      ATJUP: init_nh4sh of " << gas << " ended" << std::endl;
}
