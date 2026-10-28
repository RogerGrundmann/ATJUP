/*
 * Atmosphere General Circulation Model (ATJUP)
 * Finite difference scheme for the 3D Navier-Stokes equations on Jupiter
 * 4th-order Runge-Kutta integration
 *
 * Initial and boundary conditions for thermodynamic and chemical species fields
 */
#include "cJupiterModel.h"
#include "ATPhys.h"   // saturation_vapour_pressure
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

    // ATJUP_PDYN_INIT (2026-10-12; DEFAULT 2 on the user's word, 0 restores the zero start this
    // routine has always made): p_dyn starts in balance with the metric term of the radial
    // momentum equation, as on the other three models. rhs_u carries + (v2 + w2) / rm, the
    // centrifugal acceleration of the horizontal flow itself; with p_dyn = 0 nothing opposes it
    // and it pushes u upward in every column. The pressure solver has the term in its source,
    // but one sweep per call does not build the gradient. rhs_u takes
    // -pgrad_scale * exp_rm * dp_dyn/dr, so d p_dyn / dr = (v2 + w2) / (rm * exp_rm * pgrad_scale),
    // integrated up from zero at the bottom, trapezoid in r (inside the SeaMount v = w = 0 and
    // the pressure is carried through unchanged). The README row has the measurement.
    static const int pdyn_init = [](){ const char* e = getenv("ATJUP_PDYN_INIT"); return e ? atoi(e) : 2; }();
    // With ATJUP_METRIC_SPLIT the metric term is balanced by p_hydro, so a p_dyn that balances it
    // too would be an unopposed force: it starts at zero then (as with ATJUP_PDYN_INIT=0).
    const bool pdyn_balanced = (pdyn_init == 2) && !metric_split();
    if(pdyn_balanced){
        static const double pgrad_scale = [](){ const char* e = getenv("ATJUP_PGRAD_SCALE"); return e ? atof(e) : 1.0; }();
        double pmax = 0.0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(max:pmax)
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                double f_prev = 0.0, p_ = 0.0;
                for (int i = 0; i < im; i++) {
                    const double rm     = rad.z[i];
                    const double exp_rm = coord_stretching ? 1.0 / (rm + 1.0) : 1.0;
                    const double f = (v.x[i][j][k] * v.x[i][j][k] + w.x[i][j][k] * w.x[i][j][k]) / (rm * exp_rm * pgrad_scale);
                    if (i > 0) p_ += 0.5 * (f_prev + f) * dr;
                    p_dyn.x[i][j][k] = p_; f_prev = f;
                    if (p_ > pmax) pmax = p_;
                }
            }
        }
        printf("      ATJUP: init_PressureDynamic: p_dyn balances the metric term (v2 + w2)/rm of rhs_u; max %.6f (the model's unit) (ATJUP_PDYN_INIT=0 starts it at zero)\n", pmax);
    }

    std::cout << "      ATJUP: init_PressureDynamic ended" << std::endl;
}


// ============================================================================
// ATJUP_PDYN_GEO (on with ATJUP_GEO_BALANCE): the geostrophically balanced start, ported from
// ATNEPT's init_PressureDynamic (PDYN_GEO=2, PDYN_GEO_EXACT, PDYN_GEO_BC). Called after the
// initial boundary passes, so p_baro is built from the fields the first RK stage sees.
//
// With the radial forces carried by p_hydro (ATJUP_METRIC_SPLIT, CORIOLIS_SPLIT) what is left at
// the start is the meridional balance of rhs_v:
//     dp/dtheta = rm * Coriolis * nd_cor * 2 omega cos(theta) w + w2 cot(theta)
// (pgrad_scale 1; cot with the RK's floored sine, cos as the RK takes it). It is integrated from
// the north pole by p[j+1] = p[j-1] + 2 dthe G[j], so the centred difference rhs_v takes returns
// G exactly. G is built from the ZONAL MEAN of w over the fluid cells of each (i, j): the
// SeaMount holds w = 0, and a pressure integrated through it would differ from its neighbours in
// phi and push w along the band. p_dyn takes the barotropic part (each latitude's jet level, the
// level of the largest zonal-mean |w|); p_baro, frozen and added to p_hydro, takes what each
// level's own wind needs on top: p_baro = Q_i - p_dyn - p_hydro(start). One global constant is
// removed so that p_dyn has zero area mean.
// ============================================================================
// ATJUP_CELL_HOLD (cJupiterModel.h): the zonal means of the initial v and u over the fluid cells
// of each (i, j), frozen; the zonal Coriolis term in RHSJup acts on the departure from them.
void cJupiterModel::initCellHold(){
    if(!cell_hold()) return;
    cell_vbar0.assign((size_t)im * jm, 0.0);
    cell_ubar0.assign((size_t)im * jm, 0.0);
    double vmax = 0.0;
    for(int i = 0; i < im; i++) for(int j = 0; j < jm; j++){
        double sv = 0.0, su = 0.0; int n_ = 0;
        for(int k = 0; k < km - 1; k++) if(SeaMount.x[i][j][k] != 1.0){ sv += v.x[i][j][k]; su += u.x[i][j][k]; n_++; }
        if(n_){ cell_vbar0[(size_t)i * jm + j] = sv / n_; cell_ubar0[(size_t)i * jm + j] = su / n_; }
        vmax = std::max(vmax, std::abs(cell_vbar0[(size_t)i * jm + j]));
    }
    printf("      ATJUP: initCellHold: zonal-mean initial cells frozen, max |vbar0| %.3f m/s\n", vmax * u_0);
}


// ============================================================================
// ATJUP_VORTEX=1 (2026-10-22, default 0): the free Great Red Spot as an initial condition -- a
// geostrophically balanced elliptical anticyclone added to the balanced background (after
// initGeoBalance, so the jets' p_dyn / p_baro and the frozen cells of ATJUP_CELL_HOLD are built
// without it; the vortex's own v therefore feels the full Coriolis force). No solid body: run it with
// ATJUP_NO_SEAMOUNT=1.
//
// Pressure anomaly (barotropic up to the tropopause, then tapered to the lid -- ATJUP_VORTEX_PROFILE below)
//     p'(j, k) = P0 exp(-((k - k0)/a)^2 - ((j - j0)/b)^2)          a, b in grid cells (1 deg)
// and the wind from the model's OWN balance, term for term as RHSJup integrates it (rhs_v = rhs_w = 0
// with transport and diffusion left out), using pgrad()'s centred differences, so the start is in
// discrete balance:
//     rhs_v:  -pg (1/rm) dp/dthe + C n 2 Omega cos(the) w = 0  ->  w' = pg dp'/dthe / (rm C n 2 Omega cos)
//     rhs_w:  -pg (1/(rm s)) dp/dphi - C n 2 Omega cos(the) v = 0  ->  v' = -pg dp'/dphi / (rm s C n 2 Omega cos)
// (C = Coriolis, n = nd_cor, pg = ATJUP_PGRAD_SCALE, s and cos as the RK tables take them). Geostrophic,
// not gradient-wind: at the GRS the Rossby number is ~0.3, so the start is close to, not exactly in,
// balance. P0 is scaled so that the largest |V'| is ATJUP_VORTEX_VMAX, and its SIGN is chosen so the
// vortex is ANTICYCLONIC: in the southern hemisphere the flow on the equatorward edge goes west.
// Knobs: ATJUP_VORTEX_J (112 = 22 S), _K (180), _A (7 cells E-W half-axis ~ 8000 km at 22 S),
// _B (5 cells N-S ~ 6200 km), _VMAX (100 m/s).
// ============================================================================
void cJupiterModel::initVortex(){
    if(!vortex()) return;
    auto envd = [](const char* n, double d){ const char* e = getenv(n); return e ? atof(e) : d; };
    const int    j0 = (int)envd("ATJUP_VORTEX_J", 112), k0 = (int)envd("ATJUP_VORTEX_K", 180);
    const double a  = envd("ATJUP_VORTEX_A", 7.0), b = envd("ATJUP_VORTEX_B", 5.0);
    const double vmax = envd("ATJUP_VORTEX_VMAX", 100.0) / u_0;          // nondimensional
    const double pg = envd("ATJUP_PGRAD_SCALE", 1.0);
    static const int nd_all_ = [](){ const char* n = getenv("ATJUP_NONDIM"); return n ? atoi(n) : 0; }();
    static const int nd_cor_on_ = [](){
        const char* e = getenv("ATJUP_ND_COR"); return e ? atoi(e) : (geo_balance() ? 1 : -1); }();
    const double nd_cor = ((nd_cor_on_ >= 0 ? nd_cor_on_ : nd_all_) != 0) ? L_atm * 1.0e3 / u_0 : 1.0;
    if(nd_cor == 1.0)
        std::cout << "      ATJUP: initVortex: WARNING the Coriolis force is not at full scale (ATJUP_GEO_BALANCE=0?) -- a balanced vortex needs it" << std::endl;
    const double smin = sinthe_min();
    const double i2dthe = 1.0 / (2.0 * dthe), i2dphi = 1.0 / (2.0 * dphi);

    // shape, P0 = 1
    std::vector<double> P((size_t)jm * km, 0.0);
    for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++){
        double dk = k - k0; if(dk > km / 2) dk -= (km - 1); if(dk < -km / 2) dk += (km - 1);
        const double x = dk / a, y = (j - j0) / b;
        P[(size_t)j * km + k] = std::exp(-(x * x + y * y));
    }
    auto Pat = [&](int j, int k){ if(k < 0) k += km - 1; if(k >= km - 1) k -= km - 1; return P[(size_t)j * km + k]; };
    // wind per unit P0 at the jet level i = 0 .. im-1 (rm varies with i)
    std::vector<double> W((size_t)im * jm * km, 0.0), V((size_t)im * jm * km, 0.0);
    double sp_max = 0.0;
    for(int i = 0; i < im; i++) for(int j = 1; j < jm - 1; j++){
        const double rm = rad.z[i];
        const double s  = std::max(smin, std::abs(sin(the.z[j])));
        double c = cos(the.z[j]); if(costhe_abs() && j > 90) c = -c;
        const double f = Coriolis * nd_cor * 2.0 * omega * c;
        if(std::abs(f) < 1e-12) continue;                                 // no geostrophy on the equator
        for(int k = 0; k < km - 1; k++){
            const double dpdthe = (Pat(j + 1, k) - Pat(j - 1, k)) * i2dthe;
            const double dpdphi = (Pat(j, k + 1) - Pat(j, k - 1)) * i2dphi;
            const size_t n_ = ((size_t)i * jm + j) * km + k;
            W[n_] =  pg * dpdthe / (rm * f);
            V[n_] = -pg * dpdphi / (rm * s * f);
            sp_max = std::max(sp_max, std::sqrt(W[n_] * W[n_] + V[n_] * V[n_]));
        }
    }
    if(!(sp_max > 0.0)){ std::cout << "      ATJUP: initVortex: no wind from the pressure anomaly, skipped" << std::endl; return; }
    double P0 = vmax / sp_max;
    // anticyclonic: in BOTH hemispheres the equatorward edge flows west (w < 0) -- counterclockwise in the
    // south (the GRS), clockwise in the north; probed at the edge's steepest point, a mid level
    const int im2 = im / 2, jn = std::max(1, std::min(jm - 2, j0 - (int)std::lround(b / std::sqrt(2.0))));
    const int js = std::max(1, std::min(jm - 2, j0 + (int)std::lround(b / std::sqrt(2.0))));
    const int jeq = (j0 > 90) ? jn : js;                                  // the equatorward edge
    if(W[((size_t)im2 * jm + jeq) * km + k0] * P0 > 0.0) P0 = -P0;
    // vertical profile sz(i) of wind and p' (ATJUP_VORTEX_PROFILE, default 1; 0 = barotropic up to the lid as
    // VX2_1000): 1 up to the tropopause of the vortex row, linear to 0 at i = ATJUP_VORTEX_TOP_I (default im-3,
    // where ATJUP_BC_TOP_TAPER starts), 0 above -- the background jets' own shape above the tropopause. The
    // barotropic vortex kept its full p' in the three top layers the taper drives to rest; that unopposed
    // gradient made du/dr at i = 38 the residuum at the vortex centre (VX2_1000: 0.0135 -> 0.129 linear,
    // RES_PARTS du/dr -0.131 of the -0.131 sum). The full p' goes into p_dyn, the per-level rest
    // (sz - 1) p' into the frozen p_baro, as initGeoBalance does for the jets' per-level rest.
    // ATJUP_VORTEX_PROFILE=2: the ramp as 0.5 (1 + cos(pi x)), smooth at both ends; ATJUP_VORTEX_BASE_I moves its
    // foot below the tropopause (VTlow_300: a step 1 -> 0 made its own vertical velocity at the step -- vertical
    // viscosity on the vortex shear, the slowed wind drifting outward, Ekman pumping).
    static const int prof = [](){ const char* e = getenv("ATJUP_VORTEX_PROFILE"); return e ? atoi(e) : 1; }();
    const int itop = (int)envd("ATJUP_VORTEX_TOP_I", im - 3);
    const int tl = std::min((int)envd("ATJUP_VORTEX_BASE_I", im_tropopause[j0]), itop - 1);
    std::vector<double> sz(im, 1.0);
    if(prof) for(int i = tl + 1; i < im; i++){
        const double x = double(i - tl) / double(itop - tl);
        sz[i] = (i >= itop) ? 0.0 : (prof == 2 ? 0.5 * (1.0 + std::cos(M_PI * x)) : 1.0 - x);
    }
    // ATJUP_VORTEX_BOTTOM_I=<ib> (2026-10-24, default 0 = off: full strength down to the floor, as before): the
    // wind and p' also ramp from 0 at the floor (i = 0) to full at i = ib, same form as the upper ramp. Without it
    // the vortex stood at full strength (100 m/s) in the floor layer, where the background jets are ~0 m/s
    // (RC_3000: -4..+4 at j 112 against -45 at i 28), so its foot stayed while the upper part drifted west and
    // the columns tilted; under ATJUP_HYDRO_SPLIT p_dyn is the floor pressure and stayed with that foot. The floor
    // value sz0 = sz[0] is what goes into p_dyn, the per-level rest is (sz - sz0) p' (sz0 = 1 without the ramp,
    // so that branch is unchanged); with ATJUP_VORTEX_THERMAL the rising lower ramp is a warm anomaly under the
    // vortex, the thermal wind of a vortex strengthening upward, and the vortex's pressure lives in p_hydro.
    // ATJUP_VORTEX_LOW_I=<il> (2026-10-27, default 0 = the ramp starts at the floor, as BOTTOM_I alone): a LENS
    // vortex -- nothing (wind, p', T') below i = il, the lower ramp from 0 at il to full at ib. VB18_3000 / WK0: the
    // warm anomaly of a ramp from the floor sits in the slow deep wind (background w linear from ~0 at the floor
    // to -45 m/s at the tropopause), stays behind when the upper vortex drifts west, and regrows its own floor
    // circulation (6 -> 44 m/s) -- a fixed foot. With il in the moving layer the whole anomaly is carried by it.
    const int ib = (int)envd("ATJUP_VORTEX_BOTTOM_I", 0);
    const int il = (int)envd("ATJUP_VORTEX_LOW_I", 0);
    if(prof && ib > 0){
        if(ib >= tl) std::cout << "      ATJUP: initVortex: WARNING ATJUP_VORTEX_BOTTOM_I " << ib << " >= the base " << tl << ", lower ramp skipped" << std::endl;
        else{
            const int i0 = (il > 0 && il < ib) ? il : 0;
            if(il > 0 && il >= ib) std::cout << "      ATJUP: initVortex: WARNING ATJUP_VORTEX_LOW_I " << il << " >= ATJUP_VORTEX_BOTTOM_I " << ib << ", ignored" << std::endl;
            for(int i = 0; i < ib; i++){
                const double x = (i < i0) ? 0.0 : double(i - i0) / double(ib - i0);
                sz[i] = prof == 2 ? 0.5 * (1.0 - std::cos(M_PI * x)) : x;
            }
        }
    }
    const double sz0 = sz[0];
    // ATJUP_VORTEX_THERMAL=1 (2026-10-23, default 0; needs ATJUP_HYDRO_SPLIT and a profile): the per-level
    // rest (sz - 1) p' is carried by a TEMPERATURE anomaly instead of the frozen p_baro -- the cold top of
    // a warm-core-free anticyclone, the thermal wind of its decay with height. With the buoyancy in
    // p_hydro (ATJUP_HYDRO_SPLIT, ATJUP_BUOY_RHO_LEVEL) that pressure is the column integral of T', so it
    // is advected with the vortex instead of staying where it was put (the taper moves with it). T' from
    // the centred radial derivative of the target, (g L / u_0^2) T' / (T + T') = d((sz - 1) p')/dr per
    // cell; what the trapezoidal column integral then misses of the target goes into p_baro, so the start
    // is balanced exactly.
    static const int thermal = [](){ const char* e = getenv("ATJUP_VORTEX_THERMAL"); return e ? atoi(e) : 0; }();
    const bool th = thermal && prof && hydro_split_on() && metric_split();
    if(thermal && !th)
        std::cout << "      ATJUP: initVortex: WARNING ATJUP_VORTEX_THERMAL needs ATJUP_HYDRO_SPLIT=1, ATJUP_METRIC_SPLIT and a profile -- the rest goes into p_baro" << std::endl;
    if(prof && p_baro.empty()) p_baro.assign((size_t)im * jm * km, 0.0);
    for(int i = 0; i < im; i++) for(int j = 1; j < jm - 1; j++) for(int k = 0; k < km; k++){
        if(SeaMount.x[i][j][k] == 1.0) continue;
        const int kk = (k == km - 1) ? 0 : k;                             // km-1 is the seam copy of k = 0
        const size_t n_ = ((size_t)i * jm + j) * km + kk;
        w.x[i][j][k] += sz[i] * P0 * W[n_];
        v.x[i][j][k] += sz[i] * P0 * V[n_];
        p_dyn.x[i][j][k] += sz0 * P0 * Pat(j, kk);
        if(prof && !th) p_baro[((size_t)i * jm + j) * km + k] += (sz[i] - sz0) * P0 * Pat(j, kk);
    }
    if(th){
        const size_t n3 = (size_t)im * jm * km;
        computeBuoyancyRefLevel();
        computeHydrostaticPressure();                                     // with the vortex's wind
        std::vector<double> ph_a(n3);
        for(int i = 0; i < im; i++) for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++)
            ph_a[((size_t)i * jm + j) * km + k] = p_hydro.x[i][j][k];
        const double nd_rel = g * (L_atm * 1.0e3) / (u_0 * u_0) * envd("ATJUP_BUOY_SCALE", 1.0) * buoyancy;
        auto D = [&](int i, int j, int kk){ return (sz[std::max(0, std::min(im - 1, i))] - sz0) * P0 * Pat(j, kk); };
        double tmin = 0.0, tmax = 0.0;
        for(int i = 0; i < im; i++) for(int j = 1; j < jm - 1; j++) for(int k = 0; k < km; k++){
            if(SeaMount.x[i][j][k] == 1.0) continue;
            const int kk = (k == km - 1) ? 0 : k;
            if(Pat(j, kk) < 1e-8) continue;
            const int ia = std::max(0, i - 1), ib = std::min(im - 1, i + 1);
            const double fd = (D(ib, j, kk) - D(ia, j, kk)) / ((ib - ia) * dr);
            if(fd == 0.0) continue;
            const double T = t.x[i][j][k];
            if(!(T > 0.0) || !(buoy_ref_level[i] != 0.0)) continue;
            const double bq = g * buoy_pressure(i, j, k) / (r_mix * R_mix * T * t_ref);
            const double c = nd_rel * bq / buoy_ref_level[i];               // f' = c T' / (T + T')
            if(!(c - fd > 0.0)) continue;
            const double dT = fd * T / (c - fd);
            t.x[i][j][k] += dT;
            tmin = std::min(tmin, dT * t_ref); tmax = std::max(tmax, dT * t_ref);
        }
        computeBuoyancyRefLevel();
        computeHydrostaticPressure();
        double rmax = 0.0;
        for(int i = 0; i < im; i++) for(int j = 1; j < jm - 1; j++) for(int k = 0; k < km; k++){
            if(SeaMount.x[i][j][k] == 1.0) continue;
            const size_t q = ((size_t)i * jm + j) * km + k;
            const double r_ = D(i, j, (k == km - 1) ? 0 : k) - (p_hydro.x[i][j][k] - ph_a[q]);
            p_baro[q] += r_;
            rmax = std::max(rmax, std::abs(r_));
        }
        printf("      ATJUP: initVortex: thermal taper -- T' %+.3f .. %+.3f K, the rest the column integral misses (into p_baro) max %.3e = %.2f %% of |P0|\n",
               tmin, tmax, rmax, 100.0 * rmax / std::abs(P0));
    }
    if(prof){
        computeHydrostaticPressure();                                     // p_hydro with the vortex's p_baro
        printf("      ATJUP: initVortex: vertical profile -- full to i %d (tropopause of the row: %d), %s ramp to 0 at i %d (sz at i %d..%d:",
               tl, im_tropopause[j0], prof == 2 ? "cosine" : "linear", itop, tl, im - 1);
        for(int i = tl; i < im; i++) printf(" %.2f", sz[i]);
        printf(")\n");
        if(sz0 != 1.0){
            printf("      ATJUP: initVortex: lower ramp (ATJUP_VORTEX_BOTTOM_I %d, LOW_I %d) -- %s from 0 at i %d, sz at i 0..%d:", ib, il, prof == 2 ? "cosine" : "linear", (il > 0 && il < ib) ? il : 0, ib);
            for(int i = 0; i <= ib; i++) printf(" %.2f", sz[i]);
            printf("\n");
        }
    }
    const size_t nn = ((size_t)im2 * jm + jn) * km + k0, ns = ((size_t)im2 * jm + js) * km + k0;
    printf("      ATJUP: initVortex: anticyclone at j %d k %d (%.1f deg lat, %d E), half-axes %.1f x %.1f cells, max |V'| %.1f m/s, "
           "P0 %+.4e (%+.4f bar); w' at i %d: equatorward edge j %d %+.1f, poleward edge j %d %+.1f m/s\n",
           j0, k0, 90.0 - j0 * dthe / (M_PI / 180.0), k0, a, b, vmax * u_0, P0, P0 * p_dyn_to_bar(),
           im2, jn, P0 * W[nn] * u_0, js, P0 * W[ns] * u_0);
}


void cJupiterModel::initGeoBalance(){
    if(!pdyn_geo()) return;
    if(!metric_split()){
        std::cout << "      ATJUP: initGeoBalance: ATJUP_PDYN_GEO needs ATJUP_METRIC_SPLIT, skipped" << std::endl;
        return;
    }
    static const int nd_all_ = [](){ const char* n = getenv("ATJUP_NONDIM"); return n ? atoi(n) : 0; }();
    static const int nd_cor_on_ = [](){
        const char* e = getenv("ATJUP_ND_COR"); return e ? atoi(e) : (geo_balance() ? 1 : -1); }();
    const double nd_cor = ((nd_cor_on_ >= 0 ? nd_cor_on_ : nd_all_) != 0) ? L_atm * 1.0e3 / u_0 : 1.0;
    const double smin = sinthe_min();

    // zonal mean of w over the fluid cells, per (i, j)
    std::vector<double> wbar((size_t)im * jm, 0.0);
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 0; i < im; i++) for(int j = 0; j < jm; j++){
        double s = 0.0; int n_ = 0;
        for(int k = 0; k < km - 1; k++) if(SeaMount.x[i][j][k] != 1.0){ s += w.x[i][j][k]; n_++; }
        wbar[(size_t)i * jm + j] = n_ ? s / n_ : 0.0;
    }
    auto G = [&](int i, int j)->double{
        if(j == 0 || j == jm - 1) return 0.0;
        const double sth = std::max(smin, std::abs(sin(the.z[j])));
        double cth = cos(the.z[j]);
        if(costhe_abs() && j > 90) cth = -cth;
        const double wb = wbar[(size_t)i * jm + j];
        return rad.z[i] * Coriolis * nd_cor * 2.0 * omega * cth * wb + wb * wb * cth / sth;
    };
    auto integrate = [&](const std::vector<double>& g, std::vector<double>& P){
        P[0] = 0.0;
        for(int j = 1; j < jm; j++)
            P[j] = (j >= 2) ? P[j-2] + 2.0 * dthe * g[j-1] : P[j-1] + 0.5 * (g[j-1] + g[j]) * dthe;
    };

    // barotropic part: each latitude's jet level
    std::vector<double> g(jm), P(jm);
    std::vector<int> ia(jm, im - 1);
    for(int j = 0; j < jm; j++){
        int a = 0;
        for(int i = 1; i < im; i++) if(std::abs(wbar[(size_t)i * jm + j]) > std::abs(wbar[(size_t)a * jm + j])) a = i;
        ia[j] = a;
        g[j] = G(a, j);
    }
    integrate(g, P);
    double sw = 0.0, sp = 0.0;
    for(int j = 0; j < jm; j++){ const double s_ = std::abs(sin(the.z[j])); sw += s_; sp += s_ * P[j]; }
    const double p_off = sp / sw;
    double pmin = 0.0, pmax = 0.0;
    for(int j = 0; j < jm; j++){
        P[j] -= p_off;
        pmin = std::min(pmin, P[j]); pmax = std::max(pmax, P[j]);
        for(int i = 0; i < im; i++) for(int k = 0; k < km; k++) p_dyn.x[i][j][k] = P[j];
    }
    printf("      ATJUP: initGeoBalance: nd_cor %.1f, p_dyn = barotropic balance of the jet-level wind in rhs_v, %.6f .. %.6f (the model's unit; %.4f .. %.4f bar), jet level at the equator %d\n",
           nd_cor, pmin, pmax, pmin * p_dyn_to_bar(), pmax * p_dyn_to_bar(), ia[(jm - 1) / 2]);

    // per-level rest, frozen in p_baro
    p_baro.clear();
    computeBuoyancyRefLevel();
    computeHydrostaticPressure();                                  // p_baro still empty
    std::vector<double> pb((size_t)im * jm * km, 0.0);
    double bmin = 0.0, bmax = 0.0;
    #pragma omp parallel for schedule(static) reduction(min:bmin) reduction(max:bmax)
    for(int i = 0; i < im; i++){
        std::vector<double> gi(jm), Q(jm);
        for(int j = 0; j < jm; j++) gi[j] = G(i, j);
        integrate(gi, Q);
        for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++){
            if(SeaMount.x[i][j][k] == 1.0) continue;               // p_hydro is zero there and stays so
            const double b_ = (Q[j] - p_off) - p_dyn.x[i][j][k] - p_hydro.x[i][j][k];
            pb[((size_t)i * jm + j) * km + k] = b_;
            bmin = std::min(bmin, b_); bmax = std::max(bmax, b_);
        }
    }
    p_baro.swap(pb);
    computeHydrostaticPressure();                                  // now with p_baro
    printf("      ATJUP: initGeoBalance: p_baro (the per-level rest, added to p_hydro, frozen) %.6e .. %.6e\n", bmin, bmax);
}


// ============================================================================
// Static Pressure Initialization (hypsometric / power-law profile)
// ============================================================================
void cJupiterModel::init_PressureStatic(){
    std::cout << "\n\n\n      ATJUP: init_PressureStatic" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ATJUP_COMPOSITION: the gas constant of the mixture, R_mix, in place of R_ref = 3.75 kJ/(kg K), so that
    // p_stat is hydrostatic for the density every consumer forms, rho_local = p / (R_mix T): with R_ref the
    // two disagreed by R_ref / R_mix - 1 = 6 % (21 % against the shipped R_mix 3105). The anchor stays
    // p_ref at T = t_ref (165 K, the 1-bar level), so the deep pressures rise with the larger exponent.
    const double R_kJ = composition() ? R_mix * 1.0e-3 : R_ref;      // [kJ/(kg K)], see the unit note below
    const double exp_pressure = g / (gam * R_kJ);

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
    const double R_spec = R_kJ * 1.0e3;                // [J/(kg K)]

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

    if (composition())
        printf("      ATJUP: ATJUP_COMPOSITION: p_stat with R_mix %.1f J/(kg K) (R_ref %.1f), exponent %.4f; p at the bottom, equator %.3f bar\n",
               R_spec, R_ref * 1.0e3, exp_pressure, p_stat.x[0][(jm - 1) / 2][0]);

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
    // ATJUP_RHO_PDYN=0 (2026-10-17; default 0 with ATJUP_GEO_BALANCE, else 1): the density from
    // p_stat alone, as the buoyancy builds it (buoy_pressure). p_dyn is a potential on the constant
    // reference density r_mix; the balanced start gives it -1.25 .. +1.02 bar, which against the
    // 0.02 bar under the lid would make rho_mix negative where the chemistry reads it.
    static const double pdyn_w = [](){ const char* e = getenv("ATJUP_RHO_PDYN"); return (e ? atoi(e) != 0 : !geo_balance()) ? 1.0 : 0.0; }();
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < im; i++) {
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const double T = t.x[i][j][k] * t_ref;
                // Same guard as the buoyancy: a cell without a positive temperature has no
                // density. Written !(T > 0) so a NaN lands here instead of propagating.
                if (!(T > 0.0)) { rho_mix.x[i][j][k] = 0.0; continue; }
                // p_dyn is the nondimensional kinematic pressure, not bar — see p_dyn_to_bar()
                // in cJupiterModel.h. Added to p_stat unconverted it counted 7.79x too heavily.
                const double p_pa = (p_stat.x[i][j][k]
                                     + pdyn_w * p_dyn.x[i][j][k] * p_dyn_to_bar()) * 1.0e5;
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
    // ATJUP_BC_PHI_PERIODIC (ATPhys.h): with 1 the seam meridian k = 0 is computed like any other.
    const bool phi_per = ATPhys::phi_periodic<cJupiterModel>();
    const int  k_first = ATPhys::phi_first(phi_per);

    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 1; i < im-1; i++) {
        for (int j = 1; j < jm-1; j++) {
            const double rm       = rad.z[i];
            double sinthe         = sin(the.z[j]);
            double costhe         = cos(the.z[j]);
            if (costhe_abs() && j > 90) costhe    = -costhe;
            const double rmsinthe = rm * sinthe;

            for (int k = k_first; k < km-1; k++) {
                const int kw = ATPhys::phi_wrap(phi_per, k-1, km);
                const int ke = ATPhys::phi_wrap(phi_per, k+1, km);

                // Coriolis components
                const double Cor_r   = -2.0 * omega * sinthe * w.x[i][j][k];
                const double Cor_the = +2.0 * omega * costhe * w.x[i][j][k];
                const double Cor_phi = +2.0 * omega * (-costhe * v.x[i][j][k]
                                        + sinthe * u.x[i][j][k]);

                // Pressure-gradient components (central differences)
                const double dpdr   = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) / (2.0 * dr);
                const double dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) / (2.0 * dthe);
                const double dpdphi = (p_dyn.x[i][j][ke] - p_dyn.x[i][j][kw]) / (2.0 * dphi);

                CoriolisForce.x[i][j][k] = Coriolis * r_mix
                    * sqrt((pow(Cor_r, 2) + pow(Cor_the, 2) + pow(Cor_phi, 2)) / 3.0);

                CentrifugalForce.x[i][j][k] = centrifugal * r_mix
                    * omega * omega * rm * (1.0 + fabs(sinthe));

                // Guard the 1/t: a cell without a positive temperature has no density and
                // hence no buoyancy. bcSolidGround no longer creates 0 K cells, so this is a
                // second line of defence — one inf here used to spread through the whole
                // field within a few iterations.
                BuoyancyForce.x[i][j][k] = (t.x[i][j][k] > 0.0)
                    ? buoyancy * r_mix * g
                      * (p_stat.x[i][j][k] + p_dyn.x[i][j][k] * p_dyn_to_bar())
                      / (r_mix * R_mix * t.x[i][j][k] * t_ref) * 1e5
                    : 0.0;

                // Two conversions, both of which were missing, and this is a force DENSITY in
                // N/m3 like the three above it: p_dyn to bar (it is not bar, see p_dyn_to_bar),
                // and the nondimensional length to metres. L_atm is in KILOMETRES here, so the
                // divisor is L_atm*1e3; dividing by L_atm alone made this diagnostic 1000x too
                // large. Together the two changes leave it about 128x smaller than before.
                PresGradForce.x[i][j][k] =
                    -sqrt((pow(dpdr, 2)
                          + pow(dpdthe / rm, 2)
                          + pow(dpdphi / rmsinthe, 2)) / 3.0)
                     * p_dyn_to_bar() * 1.0e5 / (L_atm * 1.0e3);
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for Forces\n", elapsed.count() * 1e-9);
    std::cout << "      ATJUP: Forces ended" << std::endl;
}


// ============================================================================
// ATJUP_CH4_LOCAL_RHO (2026-10-26; DEFAULT 1 since the same day on the user's word, =0 restores the
// constant-density field): CH4 as a FIXED MASS FRACTION of the local density. MEASURED before the flip
// (satchk/ch4_2026-10-26/CR, PI0 restart 1000 -> 1100 vs mca_2026-10-26/M0): T, u, v, w, NH3 / H2O ice
// bit-identical -- radiation is not coupled into T by default, so only the CH4 field and the Q_rad
// diagnostic move; with ATJUP_RAD_COUPLING the methane opacity reads this field.
//
// WHY. init_vapour_cloud_ice writes every gas as magnus * r_mix * q_sat capped at r_max. Methane never
// saturates on Jupiter (t_0 90.7 K, the coldest cell ~103 K), so the cap holds everywhere and CH4 is a
// CONSTANT DENSITY, r_ch4 = 0.19 kg/m3 at the equator and 0.7 of it at the poles (137-190 g/m3 at PI0
// it 1000). Against the local density (ATJUP_LOCAL_RHO, default 1) that is a mass fraction of 0.17 at
// the bottom and above 1 over ~80 km (14 at the lid, rho 0.0135). Radiation.h reads exactly that
// q_ch4 = ch4 / rho for the methane opacity; the moist convective adjustment's density test went
// negative with it.
//
// WHAT. c = r_max(j) * rho_local / r_mix, i.e. the mass fraction the model's own reference level has
// (r_ch4 / r_mix = 0.148 at the equator, the same parabola in latitude), the density being
// computeMixtureDensity()'s; cloud and ice zero; zero inside the obstacle. Called at the end of the
// species set-up and again after a restart file is loaded, so a continued run gets the same field
// (the CH4 transport history in the file is discarded; it was the cap plus round-off). The scalar
// r_ch4 is NOT changed here: it enters r_mix and through it mue_mix / cp_mix. The abundance itself
// (0.148 by mass against 0.0142 observed) is ATJUP_COMPOSITION (setComposition below); with it set the
// mass fraction here is r_ch4 / r_mix at every latitude (methane is well mixed), without the polar 0.7.
// ============================================================================
// ============================================================================
// ATJUP_COMPOSITION=1 (2026-10-26, default 0 = off, byte-identical): the OBSERVED composition in the
// reference densities r_x that ThermalPropertiesJup turns into r_mix, c_mix, cp_mix, mue_mix, k_mix and
// R_mix, and in the deep caps of the H2O / H2S / CH4 fields.
//
// WHY. The shipped r_x are, as mass shares of r_mix = 1.2844: H2 0.673, He 0.106, H2O 0.070, CH4 0.148,
// NH3 0.0031, H2S 0.0003. r_h2 = 0.864 / r_he = 0.136 are Jupiter's VOLUME fractions entered as mass
// (He by mass is 0.230); r_ch4 = 0.19 is 10.6x the observed methane, r_h2o = 0.09 4.3x water at 2.5x
// solar, r_h2s 0.28x the observed H2S. The mixture came out with mean molar mass 2.68 (R_mix 3105)
// against 2.35, and cp_mix 10655.
//
// WHAT. Mole ratios to H2: He 0.157 (Galileo, von Zahn et al. 1998), CH4 2.37e-3 and H2S 8.9e-5
// (Galileo GPMS, Wong et al. 2004), H2O 2 x 2.5 x solar O/H = 2 x 2.5 x 4.90e-4 (Asplund et al. 2009;
// 2.5x solar = Juno MWR, Li et al. 2020), NH3 362 ppm by volume (Juno MWR deep, Li et al. 2017). Their
// mass shares w_x set r_x = w_x * r_mix with the SUM r_mix HELD at its shipped value, so p_dyn_to_bar
// and every other r_mix scaling stay put and only the composition moves; H2 takes the remainder;
// c_x = r_x / m_x. The field caps (init_vapour_cloud_ice / init_h2s) become w_x * rho_local (latitude
// parabola kept for H2O / H2S / NH3); NH3's deep plateau (t > t_add, was the constant r_nh3_add =
// 0.009 kg/m3: 0.8 % by mass at the bottom rising to 2.5 % at 2 bar, 3-10x Juno) is that same well-mixed
// share, so the NH3 mass fraction no longer rises with height below the cloud base; CH4 the uniform
// w_ch4 * rho_local (ATJUP_CH4_LOCAL_RHO follows this knob). Below the cap H2O / NH3 / H2S are
// rho_local x min(share, cold trap) (init_vapour_cloud_ice / init_h2s); the magnus * r_mix * q_sat form
// acts only with ATJUP_COMPOSITION=0 (and for CH4, which initCh4LocalRho then overwrites) -- checked
// 2026-10-27, row JUP-VAPOUR-RMIX.
// ============================================================================
// ============================================================================
// ATJUP_SVP_FIT=1 (2026-10-26, default 0 = off, byte-identical): the NH3 liquid, NH3 ice and H2O ice
// saturation vapour pressures from published curves, in the model's own four-coefficient form
//     ln E[bar] = C + (-L0 / T + del_alf ln T + del_bet T) / (1e-3 R)      (ATPhys::saturation_vapour_pressure)
// so every consumer (SaturationAdjustment, Precipitation, the moist convective adjustment, the set-up)
// takes them through the members. Not affected: the H2O liquid set (within 7 % of reference from 273 to
// 320 K), CH4, H2S, and Thermo_Jup's two-coefficient pairs, which feed only the Q_Latent diagnostic.
//
// WHY. The shipped curves are too steep: their effective latent heat L0 + del_alf T + del_bet T^2 is too
// large (NH3 ice 2550 J/g at 160 K against 1832; L0_nh3_ice was L0_nh3 x ls/lv). Against the references:
// NH3 ice x0.93 at 195.4 K, x0.18 at 160, x0.05 at 140, x0.003 at 110; NH3 liquid x1.4 at 220,
// x1.9 at the normal boiling point, x3.9 at 300; H2O ice x0.51 at 200 K, x0.09 at 150. At the NH3 ice
// deck (140-160 K) the vapour saturated 5-20x too low.
//
// WHAT (fits made 2026-10-26, satchk scratch svpfit.py):
//   H2O ice   Murphy & Koop (2005), ln p[Pa] = 9.550426 - 5723.265/T + 3.53068 ln T - 0.00728332 T --
//             the model's form exactly, no fit (611.66 Pa at 273.16 K).
//   NH3 ice   least-squares fit over 110-195.4 K to Fray & Schmitt (2009), ln p[bar] = sum_i A_i / T^i,
//             A = 15.96, -3537, -3.310e4, 1.742e6, -2.995e7; max error 0.05 %. Check: 0.0611 bar at the
//             triple point; 8.59e-4 / 2.94e-5 bar at 160 / 140 K against 8.65e-4 / 3.03e-5 from
//             Clausius-Clapeyron with 31.2 kJ/mol from the triple point.
//   NH3 liq.  least-squares fit over 195.4-330 K to the NIST Antoine equations for ammonia,
//             log10 p[bar] = 3.18757 - 506.713/(T - 80.78) below 239.6 K, 4.86886 - 1113.928/(T - 10.409)
//             above; max error 3.2 %, 5 % at 350 K. Check: 1.030 bar at 239.8 K (1.013), 10.5 bar at 300 K.
// VERIFIED 2026-10-27 (row JUP-SVP-SOURCES; the source coefficients had been written from memory): the source
// equations above are the published ones -- Murphy & Koop 2005 (QJRMS 131, 1539; the same constants in
// independent implementations), Fray & Schmitt 2009 (PSS 57, 2053; NH3-1 row of their Table, as reproduced in
// the Arrokoth ices supplement, NASA NTRS 20210019502), NIST WebBook Antoine rows for ammonia (Stull 1947,
// 164.0-239.6 and 239.6-371.5 K). The sets in ATPhys.h against them: H2O ice max 0.001 % over 110-273 K;
// NH3 ice max 0.045 % (rms 0.013 %) over 110-195.4 K; NH3 liquid max 3.2 % at 195.4 K (rms 0.9 %) over
// 195.4-330 K, 1.3 % above 239.6 K, 11.8 % at 371.5 K -- the 3.2 % is the size of NIST's own seam between its
// two rows at 239.6 K (0.993 vs 1.020 bar, 2.7 %). At the triple point the liquid set gives 0.0603 bar against
// the ice set's 0.0611 (1.3 % apart, where the two phases should agree).
// ============================================================================
void cJupiterModel::setSvpFit(){
    // The sets themselves live in the shared ATPhys.h (one copy for the four models).
    const ATPhys::SvpSet hi = ATPhys::svp_fit_h2o_ice(), ni = ATPhys::svp_fit_nh3_ice(), nl = ATPhys::svp_fit_nh3_liq();
    C_h2o_ice = hi.C; L0_h2o_ice = hi.L0; del_alf_h2o_ice = hi.del_alf; del_bet_h2o_ice = hi.del_bet;
    C_nh3_ice = ni.C; L0_nh3_ice = ni.L0; del_alf_nh3_ice = ni.del_alf; del_bet_nh3_ice = ni.del_bet;
    C_nh3     = nl.C; L0_nh3     = nl.L0; del_alf_nh3     = nl.del_alf; del_bet_nh3     = nl.del_bet;
    printf("      ATJUP: ATJUP_SVP_FIT: NH3 liquid (NIST Antoine fit), NH3 ice (Fray & Schmitt 2009 fit), H2O ice (Murphy & Koop 2005);"
           " E at 160 K: NH3 ice %.3e bar, at 200 K: H2O ice %.3e bar, at 239.8 K: NH3 liquid %.4f bar\n",
           ATPhys::saturation_vapour_pressure(160.0, C_nh3_ice, L0_nh3_ice, R_nh3, del_alf_nh3_ice, del_bet_nh3_ice),
           ATPhys::saturation_vapour_pressure(200.0, C_h2o_ice, L0_h2o_ice, R_h2o, del_alf_h2o_ice, del_bet_h2o_ice),
           ATPhys::saturation_vapour_pressure(239.8, C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3));
}


void cJupiterModel::setComposition(){
    // NH3: 362 ppm by volume (Juno MWR deep, Li et al. 2017) over an H2 mole fraction of 0.8626.
    const double x_he = 0.157, x_ch4 = 2.37e-3, x_h2s = 8.9e-5, x_h2o = 2.0 * 2.5 * 4.90e-4, x_nh3 = 362.0e-6 / 0.8626;  // per H2
    const double S = r_h2 + r_he + r_h2s + r_nh3 + r_nh4sh + r_h2o + r_ch4;
    const double w_nh4sh = r_nh4sh / S;
    // mass per mole of H2 of the mixture (NH4SH, a condensate, keeps its share)
    const double m_tot = m_h2 + x_he * m_he + x_ch4 * m_ch4 + x_h2s * m_h2s + x_h2o * m_h2o + x_nh3 * m_nh3;
    const double w_rest = 1.0 - w_nh4sh;
    const double o_h2 = r_h2, o_he = r_he, o_h2o = r_h2o, o_ch4 = r_ch4, o_h2s = r_h2s, o_nh3 = r_nh3;
    r_he  = w_rest * x_he  * m_he  / m_tot * S;
    r_ch4 = w_rest * x_ch4 * m_ch4 / m_tot * S;
    r_h2s = w_rest * x_h2s * m_h2s / m_tot * S;
    r_h2o = w_rest * x_h2o * m_h2o / m_tot * S;
    r_nh3 = w_rest * x_nh3 * m_nh3 / m_tot * S;
    r_h2  = S - (r_he + r_ch4 + r_h2s + r_h2o + r_nh3 + r_nh4sh);
    c_h2 = r_h2 / m_h2; c_he = r_he / m_he; c_ch4 = r_ch4 / m_ch4; c_h2s = r_h2s / m_h2s; c_h2o = r_h2o / m_h2o; c_nh3 = r_nh3 / m_nh3;
    printf("      ATJUP: ATJUP_COMPOSITION: mass shares of r_mix %.4f kg/m3 (held): H2 %.4f -> %.4f, He %.4f -> %.4f,"
           " H2O %.4f -> %.4f, CH4 %.4f -> %.4f, H2S %.5f -> %.5f, NH3 %.5f -> %.5f (deep plateau: this share x local"
           " density, ATJUP_R_NH3_ADD ignored)\n", S,
           o_h2 / S, r_h2 / S, o_he / S, r_he / S, o_h2o / S, r_h2o / S, o_ch4 / S, r_ch4 / S, o_h2s / S, r_h2s / S, o_nh3 / S, r_nh3 / S);
}


// ATJUP_COMPOSITION: ep_x = m_x / M_mix with M_mix = R* / R_mix, the molar mass of the gas the
// species are measured against (c / rho_local, rho_local = p / (R_mix T)), so that q_sat = ep E / p
// is a mass fraction of THAT gas. The shipped values belong to M ~ 2.22 (ep_h2o 8.125 = 18.015 / 2.217),
// neither the shipped mixture (M 2.68 from R_mix 3105) nor the corrected one (2.36): H2O and NH3
// saturation came out 6 % too high against the corrected gas. Every consumer reads the members:
// SaturationAdjustment, Precipitation, Thermo_Jup, the moist convective adjustment, the set-up.
void cJupiterModel::setCompositionEp(){
    const double R_univ = 8314.46;                      // J/(kmol K); m_x in kg/kmol, R_mix in J/(kg K)
    const double M_mix  = R_univ / R_mix;
    const double o_h2o = ep_h2o, o_nh3 = ep_nh3, o_h2s = ep_h2s, o_ch4 = ep_ch4, o_nh4sh = ep_nh4sh;
    ep_h2o = m_h2o / M_mix; ep_nh3 = m_nh3 / M_mix; ep_h2s = m_h2s / M_mix; ep_ch4 = m_ch4 / M_mix; ep_nh4sh = m_nh4sh / M_mix;
    printf("      ATJUP: ATJUP_COMPOSITION: ep = m_x / M_mix, M_mix %.4f kg/kmol: H2O %.4f -> %.4f, NH3 %.4f -> %.4f,"
           " H2S %.4f -> %.4f, CH4 %.4f -> %.4f, NH4SH %.4f -> %.4f\n", M_mix,
           o_h2o, ep_h2o, o_nh3, ep_nh3, o_h2s, ep_h2s, o_ch4, ep_ch4, o_nh4sh, ep_nh4sh);
}


void cJupiterModel::initCh4LocalRho(const char* where){
    computeMixtureDensity();
    const double d_j_half = (double)(jm - 1) / 2.0;
    double q_min = 1e30, q_max = 0.0, c_min = 1e30, c_max = 0.0;
    #pragma omp parallel for collapse(2) schedule(static) reduction(min:q_min, c_min) reduction(max:q_max, c_max)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const double par   = JupiterUtils::parabola((double)j / d_j_half);
            const double r_max_j = (R_MAX_POLE_FRAC * r_ch4 - r_ch4) * par + R_MAX_POLE_FRAC * r_ch4;
            // With ATJUP_COMPOSITION the observed abundance, well mixed: no polar reduction.
            const double frac  = composition() ? r_ch4 / r_mix : r_max_j / r_mix;
            for (int i = 0; i < im; i++) {
                const double rho = rho_mix.x[i][j][k];
                const double c   = (SeaMount.x[i][j][k] == 1.0 || !(rho > 0.0)) ? 0.0 : frac * rho;
                ch4.x[i][j][k] = c; ch4_cloud.x[i][j][k] = 0.0; ch4_ice.x[i][j][k] = 0.0;
                if (c > 0.0) { q_min = std::min(q_min, frac); q_max = std::max(q_max, frac);
                               c_min = std::min(c_min, c);    c_max = std::max(c_max, c); }
            }
        }
    }
    printf("      ATJUP: ATJUP_CH4_LOCAL_RHO (%s): CH4 = mass fraction %.4f..%.4f x local density, %.3f..%.3f g/m3; cloud and ice 0\n",
           where, q_min, q_max, c_min * 1e3, c_max * 1e3);
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
    // ATJUP_COMPOSITION (setComposition): H2O and NH3 as rho_local * min(deep share, cold trap), the
    // cold trap being the smallest saturation mass fraction ep * E(T) / p met from the bottom up to the
    // cell (ice saturation below t_0, liquid above). So the start is nowhere supersaturated -- the
    // MAGNUS_COEFF 2 x r_mix x q_sat form wrote twice saturation on the constant reference density,
    // i.e. up to ~200x saturation under the lid, which the first adjustment condensed with its latent
    // heat -- and no mass fraction rises above the tropopause. Cloud and ice start at zero; the
    // smoothing passes of the set-up are skipped for these gases (cJupiterModel.cpp), since a smoothed
    // saturated profile is supersaturated again. NH3's deep plateau r_nh3_add is replaced by the share.
    // ========================================================================
    if (composition() && (gas == "H2O" || gas == "NH3")) {
        const bool h2o_ = (gas == "H2O");
        const double Ci = h2o_ ? C_h2o_ice : C_nh3_ice, L0i = h2o_ ? L0_h2o_ice : L0_nh3_ice;
        const double ai = h2o_ ? del_alf_h2o_ice : del_alf_nh3_ice, bi = h2o_ ? del_bet_h2o_ice : del_bet_nh3_ice;
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                double q_trap = 1.0e300;
                for (int i = 0; i < im; i++) {
                    const double t_u = t.x[i][j][k] * t_ref, p_u = p_stat.x[i][j][k];
                    cloud.x[i][j][k] = 0.0; ice.x[i][j][k] = 0.0;
                    if (!(t_u > 0.0) || !(p_u > 0.0) || is_land(SeaMount, i, j, k)) { c.x[i][j][k] = 0.0; continue; }
                    const double E = (t_u >= t_0) ? ATPhys::saturation_vapour_pressure(t_u, C, L0, R, del_alf, del_bet)
                                                  : ATPhys::saturation_vapour_pressure(t_u, Ci, L0i, R, ai, bi);
                    q_trap = std::min(q_trap, ep * E / p_u);
                    const double rho_l = p_u * 1.0e5 / (R_mix * t_u);
                    c.x[i][j][k] = rho_l * std::min(r_max[j] / r_mix, q_trap);
                }
            }
        }
        auto end_c = std::chrono::high_resolution_clock::now();
        printf("      ATJUP: ATJUP_COMPOSITION: %s = rho_local x min(share %.5f, cold trap), cloud and ice 0\n",
               gas.c_str(), r_max[(jm - 1) / 2] / r_mix);
        printf(" Time measured: %.3f seconds for init_vapour_cloud_ice (%s)\n",
               std::chrono::duration_cast<std::chrono::nanoseconds>(end_c - begin).count() * 1e-9, gas.c_str());
        return;
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
                const double E_Rain = ATPhys::saturation_vapour_pressure(
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

    // ATJUP_COMPOSITION: rho_local * min(deep share, cold trap), as for H2O / NH3 in
    // init_vapour_cloud_ice; not smoothed afterwards (cJupiterModel.cpp).
    if (composition()) {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                double q_trap = 1.0e300;
                for (int i = 0; i < im; i++) {
                    const double t_u = t.x[i][j][k] * t_ref, p_u = p_stat.x[i][j][k];
                    if (!(t_u > 0.0) || !(p_u > 0.0) || is_land(SeaMount, i, j, k)) { c.x[i][j][k] = 0.0; continue; }
                    q_trap = std::min(q_trap, ep * ATPhys::saturation_vapour_pressure(t_u, C, L0, R, del_alf, del_bet) / p_u);
                    c.x[i][j][k] = p_u * 1.0e5 / (R_mix * t_u) * std::min(r_max[j] / r_mix, q_trap);
                }
            }
        }
        printf("      ATJUP: ATJUP_COMPOSITION: %s = rho_local x min(share %.5f, cold trap)\n", gas.c_str(), r_max[(jm - 1) / 2] / r_mix);
        std::cout << "      ATJUP: init_h2s of " << gas << " ended" << std::endl;
        return;
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
                const double E_Rain = ATPhys::saturation_vapour_pressure(
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
                const double E_Rain = ATPhys::saturation_vapour_pressure(
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


// ============================================================================
// ATJUP_U_INIT_CONT=<n> (2026-10-18; DEFAULT 2 since c3d0446 on the user's word, 0 = off restores; with
// ATJUP_V_MASSBAL_STRIDE 1, which keeps the closure =2 makes): THE INITIAL RADIAL
// VELOCITY FROM THE MODEL'S OWN DISCRETE CONTINUITY. Ported from ATOM_Precipitation's
// install_u_from_cells (8285348, ATM_CELLS_U_FROM_PSI), at the user's request.
//
// WHY. Since ATJUP_U_INIT_SCALE 0 (7615b7e) the initial u is exactly zero, so the prescribed
// meridional cells are horizontal sheets: v reverses with height (60 N, k = 180: +7.2 m/s at the
// bottom, -6.2 under the lid) and nothing joins the two branches. Continuity against that v asks
// for |u| up to 0.66 m/s (JNC1000, k = 180, iteration 0); the model makes 0.02 m/s by iteration
// 1000 on its own.
//
// WHAT. Per column, in the geometry the Runge-Kutta step and the pressure solver integrate in
// (metric radius rad.z, true cosine, sine floored at sinthe_min()), with D the horizontal part of
// the residuum Reporting.h prints,
//     (1/r^2) d(r^2 u)/dr = -D,     D = dv/dthe / r + cos/(r s) v + dw/dphi / (r s),
// centred in theta and phi (phi through phi_wrap), integrated upward with the trapezoid rule from
// u = 0 at the bottom wall, or at the first fluid cell over the obstacle. Solid neighbours carry
// no flux. Against the v and w the model actually holds, after the boundary passes, not against
// an analytic cell: that is ATOM's reason too.
//
// CLOSURE AT THE RIGID LID. u(lid) = 0 needs INT r^2 D dr = 0 in every column. The prescribed v
// does not have it (0.66 m/s is left at the lid at k = 180), so:
//   =1  ATOM's BOUNDED correction only: the column remainder R is taken out of the source in
//       proportion to |source|, R |S_i| / INT |S|, so no level can more than double. (ATOM's first
//       attempt scaled the two signs as its v balancing does and produced max|u| = 84 m/s next to
//       terrain.) A column more than 90 % one-signed cannot hold a closed vertical circulation and
//       keeps u = 0, as in ATOM.
//   =2  first ATOM's column balance of v (ATM_V_MASSBAL, its default since 2026-08-28):
//       v <- v - INT r v dr / INT r dr per column. The meridional part of R is the theta
//       derivative of INT r v dr, so it vanishes; what is left is the zonal part at the obstacle
//       flanks, which the correction above takes. This CHANGES v by one constant per column; the
//       shear, which is the cell, is untouched. The largest shift is printed.
//       ONLY COLUMNS THAT REACH THE BOTTOM WALL (i_topography = 0). A column standing on the
//       obstacle holds only the upper part of the cell, its mean v is that part's, and taking it
//       out REVERSED the branch: the first version, balancing every column, turned v at i = 20,
//       k = 180 from -4.90 to +9.15 m/s at 6 N and from -4.71 to +4.86 at 58 S (327 of the 444
//       sign flips in that slice were over the cone, shifts up to 23 m/s; seen by the user in the
//       zonal VTK, satchk/ucont_2026-10-18/UC2_200). Those columns keep the prescribed v and get
//       the =1 treatment.
// ATJUP_U_INIT_CONE=<m/s> (2026-10-19, DEFAULT 0.7 on the user's word, 0 restores): the one-signed columns STANDING ON THE
//       OBSTACLE (i_topography > 0) are not zeroed but closed by taking the remainder out evenly over
//       the column, src -= R / L (L its fluid height), and |u| is then clipped at the given speed.
//       WHY: =2 leaves the 5025 obstacle columns out of the v balance (it reversed v there), so
//       their source is one-signed and 3726 of them were zeroed: over the whole base footprint the
//       meridional cells had no vertical branch (i = 28, GRS2_1000: mean |u| 0.008 against 0.18 m/s
//       beside it), which the NH3 ice drew as a ring around the base. The bounded correction above
//       cannot close such a column (it removes a one-signed source entirely); the even removal can,
//       at the price of a column-uniform convergence term, and the clip keeps it from the 84 m/s
//       ATOM's unbounded first attempt made next to terrain. Clipped cells are counted.
// VOLUME, NOT MASS: ATJUP's continuity carries no density (incompressible on r_mix), where
// ATOM's is anelastic and weights with r_humid.
//
// PRINTED: max |u|, the closure left at the lid against each column's own peak |r^2 u| (ATOM's
// lesson: the endpoint, not the peak, or the check always reads ~1), the remainder removed, the
// zeroed columns, and the continuity residuum in the floored geometry before and after.
// ============================================================================
void cJupiterModel::initUFromContinuity(){
    const int mode = u_init_cont();
    std::cout << "\n      ATJUP: initUFromContinuity, ATJUP_U_INIT_CONT = " << mode << std::endl;
    const bool   periodic = ATPhys::phi_periodic<cJupiterModel>();
    const int    k_first  = ATPhys::phi_first(periodic);
    const double smin     = sinthe_min();
    auto solid = [&](int i, int j, int k){ return SeaMount.x[i][j][k] == 1.0; };
    auto s_of  = [&](int j){ return std::max(smin, std::abs(std::sin(the.z[j]))); };
    auto c_of  = [&](int j){ double c = std::cos(the.z[j]); if(costhe_abs() && j > 90) c = -c; return c; };

    // the continuity residuum in the floored geometry, max and rms over the cells the RK step
    // integrates (j = 3..jm-4) with no solid neighbour
    auto residuum = [&](double& rmax, double& rrms){
        double mx = 0.0, s2 = 0.0; long n = 0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(max:mx) reduction(+:s2,n)
        for(int i = 1; i < im-1; i++){
            for(int j = 3; j < jm-3; j++){
                const double r = rad.z[i], s = s_of(j), c = c_of(j);
                for(int k = k_first; k < km-1; k++){
                    const int kw = ATPhys::phi_wrap(periodic, k-1, km), ke = ATPhys::phi_wrap(periodic, k+1, km);
                    if(solid(i,j,k) || solid(i-1,j,k) || solid(i+1,j,k) || solid(i,j-1,k) || solid(i,j+1,k)
                       || solid(i,j,kw) || solid(i,j,ke)) continue;
                    const double res = (u.x[i+1][j][k] - u.x[i-1][j][k]) / (2.0 * dr) + 2.0 * u.x[i][j][k] / r
                        + (v.x[i][j+1][k] - v.x[i][j-1][k]) / (2.0 * dthe) / r + c / (r * s) * v.x[i][j][k]
                        + (w.x[i][j][ke] - w.x[i][j][kw]) / (2.0 * dphi) / (r * s);
                    mx = std::max(mx, std::fabs(res)); s2 += res * res; n++;
                }
            }
        }
        rmax = mx; rrms = n > 0 ? std::sqrt(s2 / n) : 0.0;
    };
    double r0max, r0rms; residuum(r0max, r0rms);

    // (=2) the column balance of v
    double vshift = 0.0;
    long nskip = 0;
    if(mode == 2) vshift = balanceColumnV(&nskip);

    // u from continuity
    double umax = 0.0;
    long ncol = 0, nzero = 0, ncone = 0, nclip = 0;
    static const double cone_lim = [](){ const char* e = getenv("ATJUP_U_INIT_CONE"); return e ? atof(e) : 0.7; }();
    std::vector<double> resid, imbal;
    #pragma omp parallel for collapse(2) schedule(static) reduction(max:umax) reduction(+:ncol,nzero,ncone,nclip)
    for(int j = 1; j < jm-1; j++){
        for(int k = k_first; k < km-1; k++){
            const int i0 = i_topography[j][k];
            for(int i = 0; i < im; i++) u.x[i][j][k] = 0.0;
            if(i0 >= im-1) continue;
            const double s = s_of(j), c = c_of(j);
            const int kw = ATPhys::phi_wrap(periodic, k-1, km), ke = ATPhys::phi_wrap(periodic, k+1, km);
            std::vector<double> src(im, 0.0);                           // r^2 D
            for(int i = i0; i < im; i++){
                const double r  = rad.z[i];
                const double vn = solid(i,j-1,k) ? 0.0 : v.x[i][j-1][k];
                const double vs = solid(i,j+1,k) ? 0.0 : v.x[i][j+1][k];
                const double ww = solid(i,j,kw)  ? 0.0 : w.x[i][j][kw];
                const double we = solid(i,j,ke)  ? 0.0 : w.x[i][j][ke];
                const double D  = (vs - vn) / (2.0 * dthe) / r + c / (r * s) * v.x[i][j][k]
                                + (we - ww) / (2.0 * dphi) / (r * s);
                src[i] = r * r * D;
            }
            double R = 0.0, S = 0.0;
            for(int i = i0; i < im-1; i++){
                R += 0.5 * (src[i] + src[i+1]) * dr;
                S += 0.5 * (std::fabs(src[i]) + std::fabs(src[i+1])) * dr;
            }
            if(!(S > 0.0)) continue;
            const double imb = std::fabs(R) / S;
            #pragma omp critical(ucont_stats)
            imbal.push_back(imb);
            const bool cone_col = cone_lim > 0.0 && i0 > 0 && imb > 0.9;
            if(imb > 0.9 && !cone_col){ nzero++; continue; }             // one-signed: no closed cell
            if(cone_col){                                                // ATJUP_U_INIT_CONE: even removal
                const double L = (im - 1 - i0) * dr;
                for(int i = i0; i < im; i++) src[i] -= R / L;
                ncone++;
            } else {
                const double corr = R / S;
                for(int i = i0; i < im; i++) src[i] -= corr * std::fabs(src[i]);
            }
            double r2u = 0.0, peak = 0.0;
            for(int i = i0; i < im-1; i++){
                r2u -= 0.5 * (src[i] + src[i+1]) * dr;
                if(i + 1 < im - 1){
                    u.x[i+1][j][k] = r2u / (rad.z[i+1] * rad.z[i+1]);
                    if(cone_col && std::fabs(u.x[i+1][j][k]) * u_0 > cone_lim){
                        u.x[i+1][j][k] = std::copysign(cone_lim / u_0, u.x[i+1][j][k]); nclip++; }
                }
                peak = std::max(peak, std::fabs(r2u));
                if(i + 1 < im - 1) umax = std::max(umax, std::fabs(u.x[i+1][j][k]));
            }
            ncol++;
            if(peak > 0.0){
                #pragma omp critical(ucont_stats)
                resid.push_back(std::fabs(r2u) / peak);                 // what the lid has to take
            }
        }
    }
    // the seam copy, and the rigid walls as bcRadius sets them
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++){
            if(periodic) u.x[i][j][km-1] = u.x[i][j][0];
        }
    for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++){ u.x[0][j][k] = 0.0; u.x[im-1][j][k] = 0.0; }

    double r1max, r1rms; residuum(r1max, r1rms);
    std::sort(resid.begin(), resid.end()); std::sort(imbal.begin(), imbal.end());
    auto q = [](const std::vector<double>& a, double f){
        return a.empty() ? 0.0 : a[std::min(a.size()-1, (size_t)(f * a.size()))]; };
    printf("      ATJUP: [U_INIT_CONT] u from continuity over %ld columns (%ld zeroed as one-signed),"
           " max |u| = %.4f m/s; v column balance %s, largest shift %.4f m/s (%ld columns on the obstacle left out)\n",
           ncol, nzero, umax * u_0, mode == 2 ? "ON" : "off", vshift * u_0, nskip);
    if(cone_lim > 0.0)
        printf("      ATJUP: [U_INIT_CONE] %ld one-signed obstacle columns closed by even removal, %ld cells clipped at %.3f m/s\n",
               ncone, nclip, cone_lim);
    printf("      ATJUP: [U_INIT_CONT] column remainder removed |R|/INT|S|: median %.4f  p95 %.4f  worst %.4f;"
           "  lid closure |r2u(lid)|/peak: median %.2e  p95 %.2e  worst %.2e\n",
           q(imbal, 0.5), q(imbal, 0.95), q(imbal, 1.0), q(resid, 0.5), q(resid, 0.95), q(resid, 1.0));
    printf("      ATJUP: [U_INIT_CONT] continuity residuum (floored geometry, nondim) max %.4e -> %.4e,"
           "  rms %.4e -> %.4e\n", r0max, r1max, r0rms, r1rms);
}


// ============================================================================
// THE COLUMN BALANCE OF v: v <- v - INT r v dr / INT r dr in every column that reaches the
// bottom wall (i_topography = 0), from j = 1 to jm-2, all k (the seam copy gets the same shift as
// k = 0 because it holds the same column). Columns on the obstacle are left out: they hold only
// the upper part of a cell, and taking out their mean reversed v there (initUFromContinuity).
// Returns the largest shift, nondimensional; *nskip counts the columns left out.
//
// Used twice: by ATJUP_U_INIT_CONT=2 at the start, and by ATJUP_V_MASSBAL_STRIDE in the loop.
// ============================================================================
//
// ATJUP_V_MASSBAL_ZONAL=1 (2026-10-22, default: follows ATJUP_VORTEX): take out only the ZONAL MEAN of
// the column means per latitude, the same shift for every column of a row. The mass constraint is
// on the latitude circle (no net meridional flux across it); a single column may carry net v,
// balanced by the zonal divergence of w -- which is all a barotropic vortex has: with the per-column
// balance VX_100 lost the whole v of its ATJUP_VORTEX anticyclone in the first 20 iterations (west /
// east edge +65 / -65 -> 0 m/s) while w kept 99.8 of 99.9 m/s. Columns on the obstacle stay out of
// the mean and unshifted, as above.
double cJupiterModel::balanceColumnV(long* nskip_out){
    static const bool zonal = [](){
        const char* e = getenv("ATJUP_V_MASSBAL_ZONAL"); return e ? atoi(e) != 0 : vortex(); }();
    double vshift = 0.0;
    long nskip = 0;
    std::vector<double> dvc(zonal ? (size_t)jm * km : 0, 0.0);
    #pragma omp parallel for collapse(2) schedule(static) reduction(max:vshift) reduction(+:nskip)
    for(int j = 1; j < jm-1; j++){
        for(int k = 0; k < km; k++){
            const int i0 = i_topography[j][k];
            if(i0 > 0){ nskip++; continue; }                             // on the obstacle: see above
            double I1 = 0.0, I0 = 0.0;
            for(int i = i0; i < im-1; i++){
                I1 += 0.5 * (rad.z[i] * v.x[i][j][k] + rad.z[i+1] * v.x[i+1][j][k]) * dr;
                I0 += 0.5 * (rad.z[i] + rad.z[i+1]) * dr;
            }
            const double dv = I1 / I0;
            if(zonal){ dvc[(size_t)j * km + k] = dv; continue; }
            for(int i = i0; i < im; i++) v.x[i][j][k] -= dv;
            vshift = std::max(vshift, std::fabs(dv));
        }
    }
    if(zonal){
        #pragma omp parallel for schedule(static) reduction(max:vshift)
        for(int j = 1; j < jm-1; j++){
            double s = 0.0; int n_ = 0;
            for(int k = 0; k < km - 1; k++) if(i_topography[j][k] == 0){ s += dvc[(size_t)j * km + k]; n_++; }
            if(!n_) continue;
            const double dv = s / n_;
            for(int k = 0; k < km; k++){
                if(i_topography[j][k] > 0) continue;
                for(int i = 0; i < im; i++) v.x[i][j][k] -= dv;
            }
            vshift = std::max(vshift, std::fabs(dv));
        }
    }
    if(nskip_out) *nskip_out = nskip;
    return vshift;
}


// ============================================================================
// ATJUP_CONE_PROJ=<n> (2026-10-18, default 0 = off and byte-identical): THE INITIAL HORIZONTAL
// WIND MADE TO FLOW AROUND THE OBSTACLE. At the user's request ("build the projection, there
// should be no flow inside the cone").
//
// WHY. The prescribed jets and cells are written into the grid and the cone is cut out of them,
// so the first fluid cells keep whatever the fields gave them: at i = 20 the wall-normal flow in
// the ring around the cone averaged 3.78 m/s against 4.76 tangential (normal/speed 0.56, the same
// at iteration 0 and 1000, REF1000 and UCM1000), up to 13.6 m/s straight into the wall. Inside
// the cone the velocity is zero; it is the cells next to it that blow into it.
//
// WHY NOT ATOM's project_initial_velocity. That runs the 3-D pressure solver with the velocity
// as source and corrects with -grad p. The shared PressureSolver.h moves information one cell per
// sweep, the cone is 40+ cells across, so the solve would need thousands of sweeps; ATOM itself
// retired its sweep knob as measured null.
//
// WHAT, per level i that cuts the obstacle. A potential-flow correction grad(phi) with
//     (1/s) d/dthe(s dphi/dthe) + (1/s^2) d2phi/dlam2 = 0      in the fluid (s floored, the RK geometry)
//     dphi/dn = -U.n                                           at every face to a solid cell
//     phi = 0                                                  at the edge of a box M cells wider than
//                                                              the obstacle's cross-section
// (M = ATJUP_CONE_PROJ_MARGIN, default 40. 10 WAS THE DEFAULT FOR ONE COMMIT (9661ce9 .. 77df90a) AND WAS TAKEN
// BACK on the user's word ("(A) back to margin 40"): at dt 0.005 over 7000 s (D5_1000) the narrow box cut the
// slowed flow next to the top -- at i = 26 a step of 17 m/s in w over one cell at the lee (west) edge, 13-23 at
// the north and south edges, while margin 40 slows w smoothly over ~20 cells -- and the continuity residuum at
// that edge rose from 0.0145 (2500 s) to 0.0181 (6700 s). The 200-it test below (280 s) was too short to show it.
// The margin test, user: "shrink the cone projection box": with the GEOM 2 base the 40-cell box reaches 7 N at
// 180 E, where the in-loop calls
// held the prescribed +40 m/s surface branch of the first Hadley cell 2-7 m/s above the same latitude
// elsewhere -- 42.5 against 35.0 at it 1000, G2_1000; 40.2 with the in-loop projection off. 200 it,
// satchk/vmax_2026-10-20: 7 N 180 E bottom 43.6 / 46.4 / 36.9 m/s for M = 40 / 20 / 10 (20 puts the box edge
// at j = 80, next to it); w off the far field in the 22 S cut beyond 14 cells from the wall 0.4-0.9 / 0.1-1.1
// / 0; wall inflow, residuum and w extrema unchanged; SOR sweeps per call ~930 -> ~220. Cost of M = 10: the
// correction is concentrated next to the body, 6-10 cells out 2.5 against 1.9 m/s, step at the box edge ~0.5.)
// on the staircase in finite-volume form: a fluid cell's faces to solid neighbours carry the
// prescribed flux and are left out of the operator, so the wall flux of the corrected field is
// zero. The correction is horizontally divergence-free and decays away from the body; it does not
// touch u. Solved by SOR per level, levels in parallel. Applied as the average of the two face
// values of the corrected velocity, so a cell next to the wall carries half of its other face.
// Then U_INIT_CONT builds u from the deflected field.
//
// IN THE LOOP (ATJUP_CONE_PROJ_STRIDE=<N>, 2026-10-19, default 10, 0 = off; in_loop = true): the same solve every N
// iterations, because nothing in the loop keeps the side walls closed (CP_200 / CP_1000: normal/speed
// back to the default 0.56 by iteration 200). Two differences from the start: phi starts from the
// previous in-loop phi of that level (warm start; the first call from zero), and the velocity is
// INCREMENTED, v += average of the two face corrections, a wall face's correction being -v (its
// face value set to zero) -- the start's form replaces the cell value by the average of its face
// values, a 1-2-1 smoothing that would act on the whole box at every call. Printed at checkpoints.
//
// WHAT IT DOES NOT DO: it turns the flow AROUND the body at each level; it does not lift any of it
// over the top (2-D per level, no radial part), and the prescribed slow-down next to the surfaces
// (bcSeaMount / bcVelSurfSur) stays. Printed: mean |U toward the wall| over the wall faces before
// and after, the rms horizontal divergence in the boxes before and after, SOR iterations.
// ============================================================================
void cJupiterModel::projectAroundObstacle(bool in_loop, bool verbose){
    if(!in_loop)
        std::cout << "\n      ATJUP: projectAroundObstacle, ATJUP_CONE_PROJ = " << cone_proj() << ", box margin " << cone_proj_margin() << " cells (ATJUP_CONE_PROJ_MARGIN)" << std::endl;
    static std::vector<std::vector<double>> phi_prev(im);              // in-loop warm start, per level
    static std::vector<std::vector<double>> phi_prev_div(im);          // the same for the ATJUP_CONE_PROJ_DIV solve
    const double smin = sinthe_min();
    std::vector<double> s(jm);
    for(int j = 0; j < jm; j++) s[j] = std::max(smin, std::abs(std::sin(the.z[j])));
    const int M = cone_proj_margin();
    // ATJUP_CONE_PROJ_DIV=<H> (2026-10-20, user: "build the local projection"; DEFAULT 3 since the same day on the
    // user's word "commit it with default 3", 0 restores): IN THE LOOP the
    // equation gets the 3D face divergence of the current flow as a source in the cells within H cells
    // (Chebyshev, 3-D) of the obstacle, so that after the correction the flux through each such cell's faces
    // balances -- with u's radial part as it is (the correction is horizontal) and zero through faces to solid
    // cells. WHY: the projection above only removes flow INTO the wall; nothing in the loop keeps continuity, and
    // over 7000 s (M40D5_1000) the flow along the base flank slowed into a staircase step with nothing taking the
    // mass: residuum at the base corner 0.0150 -> 0.0169, the same in face form (ATJUP_RES_FV, eb72849);
    // wall-adjacent rms 4.8e-3 against 1.4e-3 inside. HOW: with f = face value (average of the two cells, + the
    // correction; 0 at a wall face), s r D = sum over faces = 0 gives A phi = sum a phi_nb + s r D0, D0 the face
    // divergence of the flow before the correction (RES_FV's face form); the wall-face terms of the wall-flux
    // mode are inside D0 and drop out. Only within H cells of the body, so the correction stays harmonic and
    // decays further out (taking out the interior divergence too would put a step at the box edge). The
    // collocated grid returns the face correction to the cells as an average (the increment form below), so the
    // balance is approximate; printed: rms face divergence over those cells before and after.
    const int H = in_loop ? cone_proj_div() : 0;
    // ATJUP_CONE_PROJ_DIV: TWO SOLVES per level and call (user, 2026-10-20: "take 2."): first the divergence solve
    // above, then the wall-flux solve of the plain projection on its result (warm-started from its own previous
    // phi), both in the increment form. WHY: in face terms "cancel the wall flux" and "take out the divergence" are
    // one equation, so a single solve cannot do both on the collocated grid: the divergence solve alone left
    // 0.41 -> 0.26 m/s of wall inflow per call against 0.19 -> 0.09 without it (restart of M40D5_1000, it
    // 500 -> 600), and writing it back as the mean of the face fluxes in the near cells was worse (0.65 -> 0.58).
    // With two solves (H = 3, same restart): inflow 0.24 -> 0.063, wall-adjacent rms face divergence 1.57e-3
    // against 4.79e-3 without. Over 7000 s (PDS3_1000 against M40D5_1000, dt 0.005): residuum 0.0134 against
    // 0.0168 at 6944 s, no longer at the obstacle (the base-corner rise is gone; the max moves to 41 N 135 W,
    // 133 km); wall-adjacent rms 0.65e-3 / 0.93e-3 (centred / face) against 1.69e-3 inside; extrema equal to
    // 3-4 digits except max v at 7 N 180 E 34.9 against 38.2; ~5 % more run time.
    static std::vector<char> near_;                                     // within H cells of the body, built once
    if(H > 0 && near_.empty()){
        near_.assign((size_t)im * jm * km, 0);
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < im; i++) for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++){
            if(SeaMount.x[i][j][k] == 1.0) continue;
            bool nr = false;
            for(int di = -H; di <= H && !nr; di++){ const int ii = i + di; if(ii < 0 || ii >= im) continue;
                for(int dj = -H; dj <= H && !nr; dj++){ const int jj = j + dj; if(jj < 0 || jj >= jm) continue;
                    for(int dk = -H; dk <= H; dk++){ const int kk = k + dk; if(kk < 0 || kk >= km) continue;
                        if(SeaMount.x[ii][jj][kk] == 1.0){ nr = true; break; } } } }
            near_[((size_t)i * jm + j) * km + k] = nr ? 1 : 0;
        }
    }
    auto is_near = [&](int i, int j, int k){ return H > 0 && near_[((size_t)i * jm + j) * km + k] != 0; };
    // the 3D face divergence (RES_FV's face form) times s r, in the projection's units
    auto srD = [&](int i, int j, int k){
        auto sol = [&](int a, int b, int c){ return SeaMount.x[a][b][c] == 1.0; };
        const double r = metricRadius(rad.z[i]), rp = r + 0.5 * dr, rm = r - 0.5 * dr;
        const double sj = s[j], sp = 0.5 * (sj + s[j+1]), sm = 0.5 * (sj + s[j-1]);
        const double up = sol(i+1,j,k) ? 0.0 : 0.5 * (u.x[i][j][k] + u.x[i+1][j][k]);
        const double um = sol(i-1,j,k) ? 0.0 : 0.5 * (u.x[i][j][k] + u.x[i-1][j][k]);
        const double vp = sol(i,j+1,k) ? 0.0 : 0.5 * (v.x[i][j][k] + v.x[i][j+1][k]);
        const double vm = sol(i,j-1,k) ? 0.0 : 0.5 * (v.x[i][j][k] + v.x[i][j-1][k]);
        const double wp = sol(i,j,k+1) ? 0.0 : 0.5 * (w.x[i][j][k] + w.x[i][j][k+1]);
        const double wm = sol(i,j,k-1) ? 0.0 : 0.5 * (w.x[i][j][k] + w.x[i][j][k-1]);
        return sj * r * (rp * rp * up - rm * rm * um) / (r * r * dr) + (sp * vp - sm * vm) / dthe + (wp - wm) / dphi;
    };
    auto solid = [&](int i, int j, int k){ return SeaMount.x[i][j][k] == 1.0; };

    // wall-face statistics and box divergence, before and after
    struct Stat { double wall_sum = 0.0; long wall_n = 0; double div2 = 0.0; long div_n = 0; };
    auto stats = [&](int i, int jlo, int jhi, int klo, int khi, Stat& st){
        for(int j = jlo; j <= jhi; j++) for(int k = klo; k <= khi; k++){
            if(solid(i,j,k)) continue;
            const double vp = v.x[i][j][k], wp = w.x[i][j][k];
            if(solid(i,j-1,k)){ st.wall_sum += std::fabs(std::min(0.0,  vp)); st.wall_n++; }   // toward N: -v
            if(solid(i,j+1,k)){ st.wall_sum += std::fabs(std::max(0.0,  vp)); st.wall_n++; }   // toward S: +v
            if(solid(i,j,k-1)){ st.wall_sum += std::fabs(std::min(0.0,  wp)); st.wall_n++; }   // toward W: -w
            if(solid(i,j,k+1)){ st.wall_sum += std::fabs(std::max(0.0,  wp)); st.wall_n++; }   // toward E: +w
            if(solid(i,j-1,k) || solid(i,j+1,k) || solid(i,j,k-1) || solid(i,j,k+1)) continue;
            const double c = std::cos(the.z[j]);
            const double D = (v.x[i][j+1][k] - v.x[i][j-1][k]) / (2.0 * dthe) + c / s[j] * vp
                           + (w.x[i][j][k+1] - w.x[i][j][k-1]) / (2.0 * dphi) / s[j];
            st.div2 += D * D; st.div_n++;
        }
    };

    double t_wall0 = 0.0, t_wall1 = 0.0, t_div0 = 0.0, t_div1 = 0.0, phimax = 0.0, t_fd0 = 0.0, t_fd1 = 0.0, t_fdm = 0.0;
    long n_wall = 0, n_div = 0, n_lev = 0, it_max = 0, n_fd = 0;

    #pragma omp parallel for schedule(dynamic) reduction(+:t_wall0,t_wall1,t_div0,t_div1,n_wall,n_div,n_lev,t_fd0,t_fd1,t_fdm,n_fd) reduction(max:it_max,phimax)
    for(int i = 1; i < im-1; i++){
        int jmin = jm, jmax = -1, kmin = km, kmax = -1;
        for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++) if(solid(i,j,k)){
            jmin = std::min(jmin, j); jmax = std::max(jmax, j); kmin = std::min(kmin, k); kmax = std::max(kmax, k); }
        if(jmax < 0) continue;                                           // the level is above the body
        const int jlo = std::max(1, jmin - M), jhi = std::min(jm - 2, jmax + M);
        const int klo = std::max(1, kmin - M), khi = std::min(km - 2, kmax + M);
        const int nj = jhi - jlo + 3, nk = khi - klo + 3;                // with a halo of phi = 0
        auto at = [&](int j, int k){ return (size_t)(j - jlo + 1) * nk + (k - klo + 1); };
        Stat st0; stats(i, jlo, jhi, klo, khi, st0);

        // ATJUP_CONE_PROJ_DIV: the source, from the flow before the correction
        std::vector<double> G;
        if(H > 0){
            G.assign((size_t)nj * nk, 0.0);
            for(int j = jlo; j <= jhi; j++) for(int k = klo; k <= khi; k++){
                if(solid(i,j,k) || !is_near(i,j,k)) continue;
                const double g = srD(i, j, k);
                G[at(j,k)] = g; t_fd0 += g * g / (s[j] * s[j]); n_fd++;
            }
        }

        for(int pass = 0; pass < (H > 0 ? 2 : 1); pass++){
        const bool div_pass = H > 0 && pass == 0;                        // ATJUP_CONE_PROJ_DIV: 0 = divergence, 1 = wall flux
        std::vector<std::vector<double>>& prev = (H > 0 && pass == 0) ? phi_prev_div : phi_prev;
        std::vector<double> phi((size_t)nj * nk, 0.0);
        if(in_loop && prev[i].size() == phi.size()) phi = prev[i];

        // SOR, lexicographic within the level
        const double omega = 1.9;
        double dmax = 0.0, pmax = 0.0; int it = 0;
        for(it = 1; it <= 40000; it++){
            dmax = 0.0; pmax = 0.0;
            for(int j = jlo; j <= jhi; j++){
                const double aN = 0.5 * (s[j] + s[j-1]) / (dthe * dthe);
                const double aS = 0.5 * (s[j] + s[j+1]) / (dthe * dthe);
                const double aL = 1.0 / (s[j] * dphi * dphi);
                for(int k = klo; k <= khi; k++){
                    if(solid(i,j,k)) continue;
                    const double vp = v.x[i][j][k], wp = w.x[i][j][k];
                    double A = 0.0, B = 0.0;
                    // a fluid face enters the operator; a wall face carries the prescribed flux
                    if(div_pass && is_near(i,j,k)){
                        // ATJUP_CONE_PROJ_DIV: the wall faces are closed inside the source
                        if(!solid(i,j-1,k)){ A += aN; B += aN * phi[at(j-1,k)]; }
                        if(!solid(i,j+1,k)){ A += aS; B += aS * phi[at(j+1,k)]; }
                        if(!solid(i,j,k-1)){ A += aL; B += aL * phi[at(j,k-1)]; }
                        if(!solid(i,j,k+1)){ A += aL; B += aL * phi[at(j,k+1)]; }
                        B += G[at(j,k)];
                    } else {
                    if(solid(i,j-1,k)) B -= aN * dthe * (-vp); else { A += aN; B += aN * phi[at(j-1,k)]; }
                    if(solid(i,j+1,k)) B -= aS * dthe * ( vp); else { A += aS; B += aS * phi[at(j+1,k)]; }
                    if(solid(i,j,k-1)) B -= aL * s[j] * dphi * (-wp); else { A += aL; B += aL * phi[at(j,k-1)]; }
                    if(solid(i,j,k+1)) B -= aL * s[j] * dphi * ( wp); else { A += aL; B += aL * phi[at(j,k+1)]; }
                    }
                    if(!(A > 0.0)) continue;
                    double& p = phi[at(j,k)];
                    const double pn = p + omega * (B / A - p);
                    dmax = std::max(dmax, std::fabs(pn - p)); p = pn;
                    pmax = std::max(pmax, std::fabs(p));
                }
            }
            if(dmax <= 1.0e-10 * std::max(pmax, 1.0e-30)) break;
        }

        // the corrected velocity: the average of its two face values in each direction; a face to
        // a solid cell carries zero
        std::vector<double> vn((size_t)nj * nk, 0.0), wn((size_t)nj * nk, 0.0);
        for(int j = jlo; j <= jhi; j++){
            for(int k = klo; k <= khi; k++){
                if(solid(i,j,k)) continue;
                const double vp = v.x[i][j][k], wp = w.x[i][j][k], pp = phi[at(j,k)];
                const double fN = solid(i,j-1,k) ? 0.0 : 0.5 * (vp + v.x[i][j-1][k]) + (pp - phi[at(j-1,k)]) / dthe;
                const double fS = solid(i,j+1,k) ? 0.0 : 0.5 * (vp + v.x[i][j+1][k]) + (phi[at(j+1,k)] - pp) / dthe;
                const double fW = solid(i,j,k-1) ? 0.0 : 0.5 * (wp + w.x[i][j][k-1]) + (pp - phi[at(j,k-1)]) / (s[j] * dphi);
                const double fE = solid(i,j,k+1) ? 0.0 : 0.5 * (wp + w.x[i][j][k+1]) + (phi[at(j,k+1)] - pp) / (s[j] * dphi);
                if(!in_loop){
                    vn[at(j,k)] = 0.5 * (fN + fS);
                    wn[at(j,k)] = 0.5 * (fW + fE);
                } else {
                    const double dN = solid(i,j-1,k) ? -vp : (pp - phi[at(j-1,k)]) / dthe;
                    const double dS = solid(i,j+1,k) ? -vp : (phi[at(j+1,k)] - pp) / dthe;
                    const double dW = solid(i,j,k-1) ? -wp : (pp - phi[at(j,k-1)]) / (s[j] * dphi);
                    const double dE = solid(i,j,k+1) ? -wp : (phi[at(j,k+1)] - pp) / (s[j] * dphi);
                    vn[at(j,k)] = vp + 0.5 * (dN + dS);
                    wn[at(j,k)] = wp + 0.5 * (dW + dE);
                }
            }
        }
        if(in_loop) prev[i] = phi;
        for(int j = jlo; j <= jhi; j++) for(int k = klo; k <= khi; k++){
            if(solid(i,j,k)) continue;
            v.x[i][j][k] = vn[at(j,k)]; w.x[i][j][k] = wn[at(j,k)];
        }
        if(div_pass)                                                     // after the divergence solve alone
            for(int j = jlo; j <= jhi; j++) for(int k = klo; k <= khi; k++){
                if(solid(i,j,k) || !is_near(i,j,k)) continue;
                const double g = srD(i, j, k); t_fdm += g * g / (s[j] * s[j]);
            }
        it_max = std::max<long>(it_max, it); phimax = std::max(phimax, pmax);
        }   // pass

        Stat st1; stats(i, jlo, jhi, klo, khi, st1);
        if(H > 0)
            for(int j = jlo; j <= jhi; j++) for(int k = klo; k <= khi; k++){
                if(solid(i,j,k) || !is_near(i,j,k)) continue;
                const double g = srD(i, j, k); t_fd1 += g * g / (s[j] * s[j]);
            }
        t_wall0 += st0.wall_sum; t_wall1 += st1.wall_sum; n_wall += st0.wall_n;
        t_div0 += st0.div2; t_div1 += st1.div2; n_div += st0.div_n;
        n_lev++;
    }
    // the seam copy
    if(ATPhys::phi_periodic<cJupiterModel>())
        for(int i = 0; i < im; i++) for(int j = 0; j < jm; j++){ v.x[i][j][km-1] = v.x[i][j][0]; w.x[i][j][km-1] = w.x[i][j][0]; }

    if(!in_loop || verbose)
    printf("      ATJUP: [CONE_PROJ%s] %ld levels, SOR iterations up to %ld, max |phi| %.4e;"
           "  mean |U toward the wall| over %ld wall faces %.4f -> %.4f m/s;"
           "  rms horizontal divergence in the boxes (nondim) %.4e -> %.4e\n",
           in_loop ? "_STRIDE" : "", n_lev, it_max, phimax, n_wall, n_wall ? t_wall0 / n_wall * u_0 : 0.0, n_wall ? t_wall1 / n_wall * u_0 : 0.0,
           n_div ? std::sqrt(t_div0 / n_div) : 0.0, n_div ? std::sqrt(t_div1 / n_div) : 0.0);
    if(H > 0 && (!in_loop || verbose))
        printf("      ATJUP: [CONE_PROJ_DIV] H = %d: rms 3D face divergence x r over %ld cells next to the body %.4e"
               " -> %.4e after the divergence solve -> %.4e after the wall-flux solve\n",
               H, n_fd, n_fd ? std::sqrt(t_fd0 / n_fd) : 0.0, n_fd ? std::sqrt(t_fdm / n_fd) : 0.0, n_fd ? std::sqrt(t_fd1 / n_fd) : 0.0);
}


// ============================================================================
// ATJUP_RES_FV=1 (2026-10-20, print-only; user: test the reading of the residuum rise at the base corner):
// THE CONTINUITY RESIDUUM IN FACE FORM next to the centred form.
//
// WHY. Over 7000 s (M40D5_1000) the reported residuum rose at the bottom cells next to the west corner of
// the obstacle base (i = 1, k = 171, j = 101 / 123), all of it in dw/dphi: the flow along the flank stops in
// front of a staircase step within one cell. The reading to test: the centred differences of the steady
// query (Reporting.h) see that turn as a source although the flow, counted through the cell's faces with
// no flux into the wall, balances.
//
// WHAT. For every fluid cell i 1..im-2, j 1..jm-2, k 1..km-2: the centred residuum exactly as the steady
// query's floored line (ATJUP_RES_FLOOR; true sine if unset), and the face form
//     (r+^2 u+ - r-^2 u-) / (r^2 dr) + (s+ v+ - s- v-) / (r s dthe) + (w+ - w-) / (r s dphi)
// with each face value the average of its two cells, zero through a face to a solid cell, s+- the face
// sines as in projectAroundObstacle, r the metric radius. Printed: both at the centred maximum with the
// three directional parts, the maximum of each and where, and the rms of each over cells with a solid
// neighbour and over the rest. Changes nothing.
// ============================================================================
void cJupiterModel::resFaceDiagnostic(){
    static const double fl = [](){ const char* e = getenv("ATJUP_RES_FLOOR"); const double x = e ? atof(e) : -1.0;
                                   return (x >= 0.0 && x < 1.0) ? x : 0.0; }();
    auto solid = [&](int i, int j, int k){ return SeaMount.x[i][j][k] == 1.0; };
    auto sn = [&](int j){ return std::max(fl, std::fabs(std::sin(the.z[j]))); };
    auto cs = [&](int j){ double c = std::cos(the.z[j]); if(costhe_abs() && j > 90) c = -c; return c; };
    // the three parts of each form at (i, j, k): radial, meridional, zonal
    auto centred = [&](int i, int j, int k, double* p){
        const double r = metricRadius(rad.z[i]), s = sn(j);
        p[0] = (u.x[i+1][j][k] - u.x[i-1][j][k]) / (2.0 * dr) + 2.0 * u.x[i][j][k] / r;
        p[1] = (v.x[i][j+1][k] - v.x[i][j-1][k]) / (2.0 * dthe) / r + cs(j) / (r * s) * v.x[i][j][k];
        p[2] = (w.x[i][j][k+1] - w.x[i][j][k-1]) / (2.0 * dphi) / (r * s);
    };
    auto face = [&](int i, int j, int k, double* p){
        const double r = metricRadius(rad.z[i]), s = sn(j);
        const double rp = r + 0.5 * dr, rm = r - 0.5 * dr;
        const double up = solid(i+1,j,k) ? 0.0 : 0.5 * (u.x[i][j][k] + u.x[i+1][j][k]);
        const double um = solid(i-1,j,k) ? 0.0 : 0.5 * (u.x[i][j][k] + u.x[i-1][j][k]);
        const double vp = solid(i,j+1,k) ? 0.0 : 0.5 * (v.x[i][j][k] + v.x[i][j+1][k]);
        const double vm = solid(i,j-1,k) ? 0.0 : 0.5 * (v.x[i][j][k] + v.x[i][j-1][k]);
        const double wp = solid(i,j,k+1) ? 0.0 : 0.5 * (w.x[i][j][k] + w.x[i][j][k+1]);
        const double wm = solid(i,j,k-1) ? 0.0 : 0.5 * (w.x[i][j][k] + w.x[i][j][k-1]);
        const double sp = 0.5 * (s + sn(j+1)), sm = 0.5 * (s + sn(j-1));
        p[0] = (rp * rp * up - rm * rm * um) / (r * r * dr);
        p[1] = (sp * vp - sm * vm) / (r * s * dthe);
        p[2] = (wp - wm) / (r * s * dphi);
    };
    double cmax = 0.0, fmax = 0.0, c2w = 0.0, f2w = 0.0, c2i = 0.0, f2i = 0.0;
    long nw = 0, ni = 0;
    int ci = 0, cj = 0, ck = 0, fi = 0, fj = 0, fk = 0;
    for(int i = 1; i < im-1; i++) for(int j = 1; j < jm-1; j++) for(int k = 1; k < km-1; k++){
        if(solid(i,j,k)) continue;
        double pc[3], pf[3];
        centred(i, j, k, pc); face(i, j, k, pf);
        const double c = pc[0] + pc[1] + pc[2], f = pf[0] + pf[1] + pf[2];
        if(std::fabs(c) > cmax){ cmax = std::fabs(c); ci = i; cj = j; ck = k; }
        if(std::fabs(f) > fmax){ fmax = std::fabs(f); fi = i; fj = j; fk = k; }
        const bool wall = solid(i-1,j,k) || solid(i+1,j,k) || solid(i,j-1,k) || solid(i,j+1,k) || solid(i,j,k-1) || solid(i,j,k+1);
        if(wall){ c2w += c * c; f2w += f * f; nw++; } else { c2i += c * c; f2i += f * f; ni++; }
    }
    double pc[3], pf[3];
    centred(ci, cj, ck, pc); face(ci, cj, ck, pf);
    printf("      ATJUP: [RES_FV] centred max %.4e at i %d j %d k %d: parts radial %+.3e merid %+.3e zonal %+.3e;"
           "  face form there %+.4e: radial %+.3e merid %+.3e zonal %+.3e\n",
           cmax, ci, cj, ck, pc[0], pc[1], pc[2], pf[0] + pf[1] + pf[2], pf[0], pf[1], pf[2]);
    face(fi, fj, fk, pf); centred(fi, fj, fk, pc);
    printf("      ATJUP: [RES_FV] face-form max %.4e at i %d j %d k %d (parts %+.3e %+.3e %+.3e; centred there %+.4e);"
           "  rms centred / face: wall-adjacent cells (%ld) %.3e / %.3e, interior (%ld) %.3e / %.3e\n",
           fmax, fi, fj, fk, pf[0], pf[1], pf[2], pc[0] + pc[1] + pc[2],
           nw, nw ? std::sqrt(c2w / nw) : 0.0, nw ? std::sqrt(f2w / nw) : 0.0, ni, ni ? std::sqrt(c2i / ni) : 0.0, ni ? std::sqrt(f2i / ni) : 0.0);
}
