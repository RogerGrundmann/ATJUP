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
#include <cstdlib>   // getenv/atoi/atof for the tropopause-clamp knobs

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



// ============================================================================
// Wall-adjacent eddy viscosity around the SeaMount obstacle
// ============================================================================
//
// WHY THIS EXISTS — the measurement, so the strength is not mistaken for a free parameter.
//
// A single fluid cell beside the staircase flank of the GRS cone, (i=24, j=97, k=171) with the
// solid starting at j=98, accelerates without bound: max|w| there goes 87 -> 196 m/s over 70
// iterations, 388 m/s by iteration 320, and the run overflows shortly after. It is a purely
// local mode — the area-weighted level mean <w> and rms(w) at that level do not move (13.3 and
// 33 throughout) while the single-cell maximum triples. Controls: with the obstacle removed
// (ATJUP_NO_SEAMOUNT=1) the growth does not happen at all; replacing the radial boundary
// extrapolation (ATJUP_BC_RADIUS_COPY=1) changes nothing, bit for bit.
//
// The per-term budget of that cell (ATJUP_PROBE) names the source:
//
//   iter    w [m/s]   rhs_w   -w/(r sin) dw/dphi   -u dw/dr   diffusion   -dp/dphi
//      1      32.0     3.50           +2.44          +1.42      -0.68       -0.00
//    120     160.5    16.53          +20.26          +6.85     -10.15       -3.26
//
// The zonal gradient dw/dphi stays pinned near -20 while w grows fivefold, so -w dw/dphi is
// linear in w: an advective self-amplification with rate ~12 per nondimensional time unit,
// about 1.2 % per iteration. What should stop it is the blocking pressure a body builds up in
// front of itself, and that response is missing — dp/dphi reaches only -3.3 against +30 of
// advection. (The Poisson solver applies its one-sided obstacle stencils only when the cell
// ITSELF is solid; the fluid cell beside the wall still differences straight through the
// boundary. Fixing that is the physical repair and is NOT what this function does.)
//
// What this function does is supply the dissipation the wall region is missing. The background
// momentum diffusivity is 1/re = 1e-3, and the measured Laplacian at the probe cell is about
// -1.0e4, so 1/re contributes -10 against +30 of advection. Raising the coefficient to
// WALL_NUE_FACTOR/re over the first WALL_NUE_LAYERS cells therefore brings diffusion to the
// advection scale exactly where the flow is stagnating against a no-slip face — which is what
// an eddy viscosity does in a wall layer, and what the k-omega closure would supply here if it
// produced anything at the obstacle (nue_t is measured as 0 at this cell even with the closure
// switched on).
//
// The profile is linear in the cell distance to the nearest solid cell (Chebyshev, i.e. faces,
// edges and corners all count as distance 1), so it decays to zero at WALL_NUE_LAYERS and the
// interior solution is untouched. Both knobs are env-overridable for A/B work:
//   ATJUP_WALL_NUE=0            switches the whole treatment off (bit-identical to before)
//   ATJUP_WALL_NUE=<factor>     multiple of 1/re at the wall face itself
//   ATJUP_WALL_NUE_LAYERS=<n>   ramp depth in cells
void cJupiterModel::computeWallViscosity(){
    static const double factor = [](){
        const char* e = getenv("ATJUP_WALL_NUE");        return e ? atof(e) : 4.0; }();
    static const int    layers = [](){
        const char* e = getenv("ATJUP_WALL_NUE_LAYERS"); return e ? atoi(e) : 3;   }();

    wall_nue.initArray(im, jm, km, 0.0);
    if(factor <= 0.0 || layers <= 0){
        cout << "      ATJUP: computeWallViscosity - disabled (ATJUP_WALL_NUE=0)" << endl;
        return;
    }

    // Chebyshev distance to the nearest solid cell, capped at `layers`. Computed directly:
    // the obstacle is one compact body, so a bounded box search per cell is cheaper and far
    // simpler than a full BFS, and this runs exactly once.
    long touched = 0;
    #pragma omp parallel for collapse(2) schedule(static) reduction(+:touched)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(SeaMount.x[i][j][k] == 1.0) continue;     // inside the body: no fluid here
                int d = layers + 1;
                for(int di = -layers; di <= layers && d > 1; di++){
                    const int ii = i + di;
                    if(ii < 0 || ii >= im) continue;
                    for(int dj = -layers; dj <= layers && d > 1; dj++){
                        const int jj = j + dj;
                        if(jj < 0 || jj >= jm) continue;
                        for(int dk = -layers; dk <= layers; dk++){
                            const int kk = k + dk;
                            if(kk < 0 || kk >= km) continue;
                            if(SeaMount.x[ii][jj][kk] != 1.0) continue;
                            const int ad = std::max(std::abs(di), std::max(std::abs(dj), std::abs(dk)));
                            if(ad < d) d = ad;
                            if(d <= 1) break;
                        }
                    }
                }
                if(d > layers) continue;
                // d = 1 at the wall face -> full strength; d = layers -> just above zero.
                const double ramp = (double)(layers - d + 1) / (double)layers;
                wall_nue.x[i][j][k] = factor / re * ramp;
                touched++;
            }
        }
    }

    printf("      ATJUP: computeWallViscosity - %ld fluid cells within %d of the obstacle,"
           " peak nue_wall = %.3e (%.1f x 1/re = %.3e)\n",
           touched, layers, factor / re, factor, 1.0 / re);
}

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

    // Tropopause clamp (see the block comment below). ATJUP_IC_TROPO_CLAMP=0 restores the old
    // unbounded lapse rate for A/B; ATJUP_T_TROPOPAUSE_MIN overrides the 110 K floor, which is
    // Jupiter's observed tropopause minimum near the 0.1 bar level.
    static const bool   tropo_clamp_on     = [](){ const char* e = getenv("ATJUP_IC_TROPO_CLAMP");  return e ? atoi(e) != 0 : true;  }();
    static const double t_tropopause_min   = [](){ const char* e = getenv("ATJUP_T_TROPOPAUSE_MIN"); return e ? atof(e) : 110.0; }();

    // Local copies: im/jm/km are static const members declared without an out-of-line
    // definition, so binding them to the vector constructor's const& would odr-use them and
    // fail to link.
    const int im_v = im, jm_v = jm, km_v = km;
    i_strato_base.assign(jm_v, std::vector<int>(km_v, im_v));

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

            // Linear lapse-rate profile upward from surface, CLAMPED at the tropopause.
            //
            // The unclamped form ran the tropospheric lapse rate over the whole domain:
            // T(lid) = T(0) - gam*L_atm = T(0) - 2.0*140 = T(0) - 280 K, which put the model top
            // near 56 K against a real Jovian 110-140 K. That was the origin of the cold top —
            // measured, not inferred: the lid drifted only +1.2 K over 100 iterations, so it was
            // never a boundary-condition artifact.
            //
            // Clamping at t_tropopause_min gives an isothermal stratosphere above the level where
            // the adiabat reaches that temperature, which is the standard idealised construction
            // (and the "IC isothermal-floor" ATOM's lid pin refers to). gam is deliberately NOT
            // retuned: it also sets the pressure exponent g/(gam*R_ref) used here and in
            // SaturationAdjustmentJup, so changing it would move the whole p_stat profile and the
            // cloud bases with it.
            //
            // The clamp height follows latitude for free: z_c = (T(0) - t_tropopause_min)/gam, and
            // T(0) is largest at the equator, so the tropopause sits higher there. That reproduces
            // the sense of initTropopauseLayers() (125 km equator / 115 km pole) without coupling
            // to it — with the default 110 K the crossing lands near 113 km.
            const int i_clamp_none = im;
            int i_clamp = i_clamp_none;
            for (int i = 1; i < im; i++) {
                const double height = get_layer_height(i);
                const double t_lapse = -gam * height + t.x[0][j][k] * t_ref;   // [K]
                if (tropo_clamp_on && t_lapse < t_tropopause_min) {
                    t.x[i][j][k] = t_tropopause_min / t_ref;
                    if (i_clamp == i_clamp_none) i_clamp = i;
                } else {
                    t.x[i][j][k] = t_lapse / t_ref;
                }
            }
            i_strato_base[j][k] = i_clamp;
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

    // p_stat ~ (T/t_ref)^(g/(gam*R_ref)) is the POLYTROPIC relation: it assumes the lapse rate is
    // exactly gam. That holds through the troposphere, but init_temperature now clamps the profile
    // to an isothermal stratosphere above i_strato_base, and there the polytropic form would hold
    // p_stat CONSTANT with height — pressure would stop decreasing, which is unphysical and would
    // corrupt every p_stat consumer (CIA opacity ~ P^2/T, cloud bases, saturation ratios).
    //
    // So above the clamp switch to the isothermal hydrostatic law,
    //   p(z) = p_c * exp(-(z - z_c)/H),   H = R_spec*T_iso/g,
    // anchored on the last adiabatic layer. R_ref is in kJ/(kg K) (g/(gam*R_ref) = 25.92/(2*3.75)
    // = 3.456 only works out with gam in K/km and R_ref in kJ/(kg K)), so the specific gas
    // constant is R_ref*1e3 = 3750 J/(kg K); with T_iso = 110 K that gives H ~ 15.9 km.
    //
    // Below the clamp the analytic formula is kept untouched, so the troposphere — and every
    // calibration already tuned against it — is bit-identical.
    const double R_spec = R_ref * 1.0e3;               // [J/(kg K)]

    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const int isb = (i_strato_base.empty() ? im : i_strato_base[j][k]);

            const int i_adiabatic_end = (isb < im) ? isb : im;
            for (int i = 0; i < i_adiabatic_end; i++) {
                const double t_u = t.x[i][j][k] * t_ref;
                p_stat.x[i][j][k] = p_ref * pow(t_u / t_ref, exp_pressure);
            }

            if (isb < im && isb > 0) {
                const double t_iso = t.x[isb][j][k] * t_ref;                  // [K]
                const double H     = R_spec * t_iso / g;                      // [m]
                const double p_c   = p_stat.x[isb-1][j][k];                   // anchor
                const double z_c   = (double)get_layer_height(isb-1) * 1.0e3; // [m]
                for (int i = isb; i < im; i++) {
                    const double z = (double)get_layer_height(i) * 1.0e3;     // [m]
                    p_stat.x[i][j][k] = p_c * exp(-(z - z_c) / H);
                }
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
// Local density of the gas mixture, rho = p/(R_mix*T) — the ideal gas law, with p converted
// from bar to Pa exactly as the buoyancy expression in Forces() does it.
//
// It was NEVER COMPUTED. resetArrays initialises rho_mix twice, to 0.0 and then to 1.0, and no
// line anywhere in the model assigns to it again — the only occurrence of rho_mix.x[][][] in the
// whole codebase is a READ, in the NH4SH Stokes settling velocity in RHS_Jup_Turb.cpp. So the
// field sat at exactly 1.0 kg/m3 in every cell for the entire run, which is what a zonal slice
// at iteration 200 showed: min = max = mean = 1.000000 over all 7421 cells.
//
// The consequence for the physics was small, which is why it survived: the settling velocity
// goes as (rg_nh4sh - rho_mix) with rg_nh4sh = 1170 kg/m3, so using 1.0 instead of the true
// 0.1..1.3 kg/m3 was at most a 0.1 % error there. The consequence for anyone LOOKING at the
// model was not small — rho_mix is written to all three ParaView slices, where it was a flat
// constant masquerading as a computed field.
//
// Placed next to Forces() and called with it, so it shares the physics-block cadence and is
// consistent with the p_stat and t it is built from. Whole grid, boundaries included, unlike
// Forces() itself, since nothing here needs neighbours.
void cJupiterModel::computeMixtureDensity(){
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < im; i++) {
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const double T = t.x[i][j][k] * t_ref;
                // Same guard as the buoyancy: a cell without a positive temperature has no
                // density. Written !(T > 0) so a NaN lands here instead of propagating.
                if (!(T > 0.0)) { rho_mix.x[i][j][k] = 0.0; continue; }
                const double p_pa = (p_stat.x[i][j][k] + p_dyn.x[i][j][k]) * 1.0e5;
                const double rho  = p_pa / (R_mix * T);
                rho_mix.x[i][j][k] = std::isfinite(rho) ? rho : 0.0;
            }
        }
    }
}

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

                // Guard the 1/t: a cell without a positive temperature has no density and
                // hence no buoyancy. bcSolidGround no longer creates 0 K cells, so this is a
                // second line of defence — one inf here used to spread through the whole
                // field within a few iterations.
                BuoyancyForce.x[i][j][k] = (t.x[i][j][k] > 0.0)
                    ? buoyancy * r_mix * g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])
                      / (r_mix * R_mix * t.x[i][j][k] * t_ref) * 1e5
                    : 0.0;

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

    if (gas == "CH4") {
        r_max_equator      = r_ch4;
        r_max_pole         = R_MAX_POLE_FRAC * r_ch4;
        cloud_loc_equator  = 1.0;
        cloud_loc_pole     = 0.0;
    }
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
