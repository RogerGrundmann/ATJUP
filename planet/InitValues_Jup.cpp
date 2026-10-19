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
    long ncol = 0, nzero = 0;
    std::vector<double> resid, imbal;
    #pragma omp parallel for collapse(2) schedule(static) reduction(max:umax) reduction(+:ncol,nzero)
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
            if(imb > 0.9){ nzero++; continue; }                          // one-signed: no closed cell
            const double corr = R / S;
            for(int i = i0; i < im; i++) src[i] -= corr * std::fabs(src[i]);
            double r2u = 0.0, peak = 0.0;
            for(int i = i0; i < im-1; i++){
                r2u -= 0.5 * (src[i] + src[i+1]) * dr;
                if(i + 1 < im - 1) u.x[i+1][j][k] = r2u / (rad.z[i+1] * rad.z[i+1]);
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
double cJupiterModel::balanceColumnV(long* nskip_out){
    double vshift = 0.0;
    long nskip = 0;
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
            for(int i = i0; i < im; i++) v.x[i][j][k] -= dv;
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
//     phi = 0                                                  at the edge of a box 40 cells wider than
//                                                              the obstacle's cross-section
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
        std::cout << "\n      ATJUP: projectAroundObstacle, ATJUP_CONE_PROJ = " << cone_proj() << std::endl;
    static std::vector<std::vector<double>> phi_prev(im);              // in-loop warm start, per level
    const double smin = sinthe_min();
    std::vector<double> s(jm);
    for(int j = 0; j < jm; j++) s[j] = std::max(smin, std::abs(std::sin(the.z[j])));
    const int M = 40;
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

    double t_wall0 = 0.0, t_wall1 = 0.0, t_div0 = 0.0, t_div1 = 0.0, phimax = 0.0;
    long n_wall = 0, n_div = 0, n_lev = 0, it_max = 0;

    #pragma omp parallel for schedule(dynamic) reduction(+:t_wall0,t_wall1,t_div0,t_div1,n_wall,n_div,n_lev) reduction(max:it_max,phimax)
    for(int i = 1; i < im-1; i++){
        int jmin = jm, jmax = -1, kmin = km, kmax = -1;
        for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++) if(solid(i,j,k)){
            jmin = std::min(jmin, j); jmax = std::max(jmax, j); kmin = std::min(kmin, k); kmax = std::max(kmax, k); }
        if(jmax < 0) continue;                                           // the level is above the body
        const int jlo = std::max(1, jmin - M), jhi = std::min(jm - 2, jmax + M);
        const int klo = std::max(1, kmin - M), khi = std::min(km - 2, kmax + M);
        const int nj = jhi - jlo + 3, nk = khi - klo + 3;                // with a halo of phi = 0
        auto at = [&](int j, int k){ return (size_t)(j - jlo + 1) * nk + (k - klo + 1); };
        std::vector<double> phi((size_t)nj * nk, 0.0);
        if(in_loop && phi_prev[i].size() == phi.size()) phi = phi_prev[i];

        Stat st0; stats(i, jlo, jhi, klo, khi, st0);

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
                    if(solid(i,j-1,k)) B -= aN * dthe * (-vp); else { A += aN; B += aN * phi[at(j-1,k)]; }
                    if(solid(i,j+1,k)) B -= aS * dthe * ( vp); else { A += aS; B += aS * phi[at(j+1,k)]; }
                    if(solid(i,j,k-1)) B -= aL * s[j] * dphi * (-wp); else { A += aL; B += aL * phi[at(j,k-1)]; }
                    if(solid(i,j,k+1)) B -= aL * s[j] * dphi * ( wp); else { A += aL; B += aL * phi[at(j,k+1)]; }
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
        if(in_loop) phi_prev[i] = phi;
        for(int j = jlo; j <= jhi; j++) for(int k = klo; k <= khi; k++){
            if(solid(i,j,k)) continue;
            v.x[i][j][k] = vn[at(j,k)]; w.x[i][j][k] = wn[at(j,k)];
        }

        Stat st1; stats(i, jlo, jhi, klo, khi, st1);
        t_wall0 += st0.wall_sum; t_wall1 += st1.wall_sum; n_wall += st0.wall_n;
        t_div0 += st0.div2; t_div1 += st1.div2; n_div += st0.div_n;
        n_lev++; it_max = std::max<long>(it_max, it); phimax = std::max(phimax, pmax);
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
}
