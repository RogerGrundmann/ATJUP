/*
 * Atmosphere General Circulation Modell(ATJUP) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in aa spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the boundary and initial conditions for diverse variables
*/
#include "cJupiterModel.h"
#include "Utils.h"
#include "ATPhys.h"   // clausius_clapeyron(): the saturation formula that
                                       // matches ATJUP's coeff_*_A/B pairs

using namespace std;
using namespace JupiterUtils;

void cJupiterModel::TropopauseLocation(){
//    cout << endl << "      ATJUP: TropopauseLocation" << endl;

// parabolic tropopause location distribution from pole to pole assumed
//
// WHERE THE END POINTS COME FROM (ATJUP_TROPO_FROM_CONFIG, default 1).
//
// im_tropopause is what bounds the overturning cells vertically: VelocityInitializerJup::init_u
// reads nothing else to decide how high the radial branch of the Hadley/Ferrel cells reaches,
// and init_v_or_w uses it to split the linear profile from the above-tropopause decay. The
// hard-coded end points are i_max = 40 and i_beg = 30, labelled "about 125 km" and "about 100 km"
// in cJupiterModel.h -- but a layer on this grid is L_atm/(im-1) = 140/40 = 3.5 km, so 40 is
// 140 km, the MODEL LID, and 30 is 105 km. Measured in the initial state at j=112, k=90: the
// ramp peaks at 19.4 m/s at 91 km and returns to zero only at i=40, i.e. one triangular branch
// spanning the entire shell.
//
// The heights are already configured, in km, as tropopause_equator (125.0) and tropopause_pole
// (115.0). Those had no effect on this: init_tropopause_layers() computes tropopause_layers from
// them, but the getter at cJupiterModel.h:118 overwrites it with im_tropopause[j] before anyone
// reads it, so the hard-coded pair won. They now set the end points, looked up against the actual
// layer heights rather than assuming a spacing, so the construction survives coord_stretching.
// On the default grid this gives i_max = 36 (126 km) and i_beg = 33 (115.5 km).
//
// This bounds the branch below the lid; it does NOT make the cells shorter than the troposphere.
// init_u has no vertical stacking at all -- one ramp per latitude, up to 2/3 of the tropopause
// height and back down -- and the cell structure comes entirely from the sign alternation across
// latitude in u_amplitude. A circulation that closes lower than the tropopause needs a different
// ramp, which is a decision about the intended physics.
//
// ATJUP_TROPO_FROM_CONFIG=0 restores the hard-coded 40/30 for A/B work.
    static const int from_config = [](){
        const char* e = getenv("ATJUP_TROPO_FROM_CONFIG"); return e ? atoi(e) : 1; }();

    auto layer_for_height = [&](double h_km){
        int best = 0; double best_d = 1.0e30;
        for(int i = 0; i < im; i++){
            const double d = std::fabs((double)get_layer_height(i) - h_km);
            if(d < best_d){ best_d = d; best = i; }
        }
        return best;
    };

    const int i_max_eff = from_config ? layer_for_height(tropopause_equator) : i_max;
    const int i_beg_eff = from_config ? layer_for_height(tropopause_pole)    : i_beg;

    printf("      ATJUP: tropopause equator = layer %d (%.1f km), pole = layer %d (%.1f km)%s\n",
           i_max_eff, (double)get_layer_height(i_max_eff),
           i_beg_eff, (double)get_layer_height(i_beg_eff),
           from_config ? "" : "  [hard-coded 40/30, ATJUP_TROPO_FROM_CONFIG=0]");

    im_tropopause = std::vector<int>(jm, 0);
    int j_half = (jm-1)/2;
    double d_j_half = (double)j_half;
    double trop_u2_eff = (double)(i_beg_eff - i_max_eff);
//    double trop_u2_eff = (double)(i_beg_trop - i_max_trop);
// computation of the tropopause from pole to pole
    #pragma omp parallel for
    for(int j = 0; j < jm; j++){
        double d_j = (double)j;
        im_tropopause[j] = (int)((trop_u2_eff * (d_j * d_j/(d_j_half * d_j_half)
            - 2.0 * d_j/d_j_half)) + (double)i_beg_eff);
//            - 2.0 * d_j/d_j_half)) + (double)i_beg_trop);
// cout << "   j = " << j << "   im_tropopause[j] = " << im_tropopause[j] << endl;
    }

//    cout << "      ATJUP: TropopauseLocation ended" << endl;
    return;
}
/*
*
*/
void cJupiterModel::Latent_Heat(){
    cout << endl << "      ATJUP: Latent_Heat" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Q_Latent.x[0][j][k] = 0.0;
            Q_Sensible.x[0][j][k] = 0.0;  // sensible heat in [W/m2] from energy transport equation
        }
    }

    #pragma omp parallel for
    for(int j = 1; j < jm-1; j++){
        double sinthe = sin(the.z[j]);
        for(int k = 1; k < km-1; k++){
            for(int i = im-2; i >= 1; i--){
                double rm = rad.z[i];
                double rmsinthe = rm * sinthe;
                double exp_rm = coord_stretching ? 1.0 / (rm + 1.0) : 1.0;
                double Latency_Ice = 0.0;

                double t_u = t.x[i][j][k] * t_ref;
                double p_u = p_stat.x[i][j][k];

                // Saturation vapour pressures [bar], Clausius-Clapeyron: E = exp(A/T + B).
                //
                // These used to call JupiterUtils::exp_func, which is ATOM's MAGNUS/TETENS form
                //     exp(c1 * (T - 273.15) / (T - c2))
                // and expects ATOM's terrestrial coefficient pair (c1 ~ 17.27, c2 ~ 35.86 K).
                // ATJUP's coeff_*_A/B are Clausius-Clapeyron coefficients instead (A = -L/R_v:
                // -4961.04 K for H2O, i.e. L = 2.29e6 J/kg; -2836.56 K for NH3), which is the
                // pair ATPhys::clausius_clapeyron() consumes. Feeding them to
                // the Magnus form gave, at the Jovian t_u ~ 110 K of the upper layers,
                //     -4961.04 * (110 - 273.15) / (110 - 13.07) = +8350   ->   exp(8350) = inf,
                // and the mixing ratio below then evaluated inf/(p - inf) = inf/-inf = NaN.
                // That was the FIRST invalid floating-point operation in the whole run (caught
                // with feenableexcept(FE_INVALID) under gdb) and it poisoned Q_Latent.
                //
                // Units: exp(A/T + B) is in BAR (check: H2O at 273.15 K gives 6.1e-3 bar =
                // 6.1 hPa, the textbook value; NH3 at its 195.5 K triple point gives 0.062 bar
                // against the tabulated 0.0606 bar), and p_stat is in bar, so the 1e3 factor
                // that used to convert to hPa has to go with it — the ratio E/p must be formed
                // in one consistent unit.
                double E_Rain = ATPhys::clausius_clapeyron(t_u, coeff_h2o_A, coeff_h2o_B);
                double E_Ice = ATPhys::clausius_clapeyron(t_u, coeff_h2o_A_i, coeff_h2o_B_i);
                // Saturation mixing ratio q = ep*E/(p - E). Guard the denominator: above the
                // critical point, and in the solid SeaMount cells where p_stat is 0, p - E can
                // vanish or go negative, which is not a physical state but must not produce
                // inf/NaN. Clamp to a small positive residual pressure.
                double q_Rain = ep_h2o * E_Rain / std::max(p_u - E_Rain, 1.0e-12);  // h2o vapour amount at saturation with water formation in kg/kg
                double q_Ice = ep_h2o * E_Ice / std::max(p_u - E_Ice, 1.0e-12);  // h2o vapour amount at saturation with ice formation in kg/kg

                double E_Rain_nh3 = ATPhys::clausius_clapeyron(t_u, coeff_nh3_A, coeff_nh3_B);
                double E_Ice_nh3 = ATPhys::clausius_clapeyron(t_u, coeff_nh3_A_i, coeff_nh3_B_i);
                double q_Rain_nh3 = ep_nh3 * E_Rain_nh3 / std::max(p_u - E_Rain_nh3, 1.0e-12);  // nh3 vapour amount at saturation with water formation in kg/kg
                double q_Ice_nh3 = ep_nh3 * E_Ice_nh3 / std::max(p_u - E_Ice_nh3, 1.0e-12);  // nh3 vapour amount at saturation with ice formation in kg/kg

                double u_av = 0.5 * (u.x[i+1][j][k] + u.x[i-1][j][k]);
                double v_av = 0.5 * (v.x[i+1][j][k] + v.x[i-1][j][k]);
                double w_av = 0.5 * (w.x[i+1][j][k] + w.x[i-1][j][k]);

                double velocity_av =
                        - sqrt((pow(u_av, 2)
                        + pow(v_av, 2)
                        + pow(w_av, 2))/3.0);

                double dtdr = (t.x[i+1][j][k] - t.x[i-1][j][k])/(2.0 * dr) * exp_rm;
                double dtdthe = (t.x[i][j+1][k] - t.x[i][j-1][k])/(2.0 * rm * dthe);
                double dtdphi = (t.x[i][j][k+1] - t.x[i][j][k-1])/(2.0 * rmsinthe * dphi);

                double dh2odr = (h2o.x[i+1][j][k] - h2o.x[i-1][j][k])/(2.0 * dr) * exp_rm;
                double dh2odthe = (h2o.x[i][j+1][k] - h2o.x[i][j-1][k])/(2.0 * rm * dthe);
                double dh2odphi = (h2o.x[i][j][k+1] - h2o.x[i][j][k-1])/(2.0 * rmsinthe * dphi);

                double dnh3dr = (nh3.x[i+1][j][k] - nh3.x[i-1][j][k])/(2.0 * dr) * exp_rm;
                double dnh3dthe = (nh3.x[i][j+1][k] - nh3.x[i][j-1][k])/(2.0 * rm * dthe);
                double dnh3dphi = (nh3.x[i][j][k+1] - nh3.x[i][j][k-1])/(2.0 * rmsinthe * dphi);

                double dtemp = dtdr + dtdthe + dtdphi;
                double dh2o = dh2odr + dh2odthe + dh2odphi;
                double dnh3 = dnh3dr + dnh3dthe + dnh3dphi;


                if(h2o.x[i][j][k] >= q_Rain)  
                    Q_Latent.x[i][j][k] = lv_h2o * velocity_av * dh2o * u_0/(L_atm * (im-1));
                else  Q_Latent.x[i][j][k] = 0.0;

                if(h2o.x[i][j][k] >= q_Ice)  
                    Latency_Ice = ls_h2o * velocity_av * dh2o * u_0/(L_atm * (im-1));
                else  Latency_Ice = 0.0;



                // The NH3 contribution ACCUMULATES onto the H2O one. These two branches used to
                // read `else Q_Latent = 0.0` / `else Latency_Ice = 0.0`, which threw the H2O
                // term away wherever NH3 happened to be subsaturated — and in the deep, warm
                // layers NH3 always is (its saturation pressure at 325 K is ~20 bar against an
                // ambient ~11 bar, so it cannot condense there at all). The result was a
                // Q_Latent that was zero over most of the domain. A species that does not
                // condense contributes nothing; it does not erase the species that does.
                if(nh3.x[i][j][k] >= q_Rain_nh3)
                    Q_Latent.x[i][j][k] = Q_Latent.x[i][j][k] + lv_nh3
                        * velocity_av * dnh3 * u_0/(L_atm * (im-1));

                if(nh3.x[i][j][k] >= q_Ice_nh3)
                    Latency_Ice = Latency_Ice + ls_nh3 * velocity_av * dnh3 * u_0/(L_atm * (im-1));



                Q_Latent.x[i][j][k] = Q_Latent.x[i][j][k] + Latency_Ice;  // latent heat in [W/m³] from energy transport equation



                Q_Sensible.x[i][j][k] = r_mix * cp_mix 
                    * velocity_av * u_0 * dtemp * t_ref/pow(L_atm * (im-1),2);  // sensible heat in [W/m³] from energy transport equation


                if(SeaMount.x[i][j][k] == 1.0){
                    Q_Latent.x[i][j][k] = 0.0;
                    Q_Sensible.x[i][j][k] = 0.0;
                }
            }
        }
    }


    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for Latent_Heat\n", elapsed.count() * 1e-9);

    cout << "      ATJUP: Latent_Heat ended" << endl;
    return;
}
/*
*
*/
