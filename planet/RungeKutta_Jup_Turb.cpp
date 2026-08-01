/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * 4th order Runge-Kutta integration of the prognostic fields assembled in RHS_Jup_Turb.cpp
 *
 * Renamed from RungeKutta_Jup.cpp when the turbulent kinetic energy k* and its dissipation
 * dis* (epsilon* for k-epsilon, omega* for k-omega / SST) joined the integrated set, mirroring
 * ATOM_Precipitation/atmosphere/RungeKutta_Atm_Turb.cpp. Previously TurbulenceJup::advance()
 * integrated them on its own with a Patankar splitting and without any transport term; they are
 * now two more scalars of this RK4 system, so they see the same advection, the same turbulent
 * diffusion and the same four-stage time integration as temperature and the species.
*/

#include "cJupiterModel.h"

#include <cstdint>
#include <cstring>

using namespace std;

// Boussinesq base state for the buoyancy term.
//
// ATJUP's rhs_u carried the ABSOLUTE buoyancy g*(p_stat+p_dyn)/(r_mix*R_mix*t*t_ref). That
// quantity is positive in every cell and of order g, so the momentum equation received a
// systematic upward acceleration everywhere, opposed only by the radial pressure gradient.
// Any residual between the two accumulates, and it did: the vertical velocity grew monotonically
// by roughly 1 m/s per iteration (max w = 130 m/s at iter 60, 209 at 150, 309 at 240, 423 at
// 376) until the deepest three levels overflowed to NaN at iteration 366 and the whole field
// collapsed. Values of several hundred m/s are unphysical for a Jovian vertical wind long
// before that.
//
// The fix is ATOM's (cAtmosphereModel::computeLevelMeanTemperature + the (t - t_ref_level[i])
// anomaly in RHS_Atm_Turb.cpp): subtract the horizontal mean at each height, so the body force
// has ZERO mean at every level and only horizontal density contrasts drive vertical motion.
// The mean is what hydrostatic balance carries, and it is not the flow's job to fight it.
//
// Unlike ATOM this keeps ATJUP's OWN expression and units rather than importing ATOM's
// empirical g*dt/u_0 coefficient — the anomaly is formed from exactly the term that was there
// before, so nothing needs recalibrating; only its horizontal mean is removed.
//
// The mean is area-weighted with sin(theta) (spherical area element) and taken over fluid
// cells only; a level that is entirely inside the SeaMount gets 0, which leaves those cells
// with their own value and is harmless since RungeKuttaJup skips them anyway.
void cJupiterModel::computeBuoyancyRefLevel(){
    if((int)buoy_ref_level.size() != im) buoy_ref_level.assign(im, 0.0);

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < im; i++){
        double sum = 0.0, wsum = 0.0;
        for(int j = 0; j < jm; j++){
            const double wgt = sin(the.z[j]);              // spherical area weight
            for(int k = 0; k < km; k++){
                if(SeaMount.x[i][j][k] == 1.0) continue;   // solid cell: no fluid state
                if(!(t.x[i][j][k] > 0.0)) continue;        // also catches NaN
                // Same pressure the buoyancy itself uses — the mean has to be the mean OF the
                // quantity whose anomaly is taken, or the anomaly no longer has zero mean.
                const double b = g * buoy_pressure(i, j, k)
                               / (r_mix * R_mix * t.x[i][j][k] * t_ref);
                if(!std::isfinite(b)) continue;
                sum  += wgt * b;
                wsum += wgt;
            }
        }
        buoy_ref_level[i] = (wsum > 0.0) ? sum / wsum : 0.0;
    }
}

// ---------------------------------------------------------------------------------------------
// Hydrostatic pressure perturbation: split the buoyancy so the pressure equation is asked only
// for what it can actually deliver.
//
// THE PROBLEM THIS SOLVES. With the body forces in their proper units the buoyancy is 8 to 13 in
// the units rhs_u is written in, against transport terms of order one, and it points radially. In
// a real atmosphere almost all of that is carried by the hydrostatic pressure gradient and never
// accelerates anything; only the small residual drives vertical motion. This model expected
// PressureSolverJup to produce that balance from the projection. Measured on a restart, it does
// not: the pressure cancels 55-63 % of the buoyancy low down, nothing by 80 km, and above 90 km it
// ADDS to it. What is left over, 4 to 8, is what accelerates the radial wind to 300 m/s in 250
// iterations until the run dies.
//
// It is not a solver-convergence problem — max|u| at iteration 60 is 139 with one relaxation
// sweep, 93 with fifty and 90 with two hundred, so it has saturated. The projection cannot build a
// hydrostatic pressure because p_dyn is given a zero normal gradient at BOTH radial walls, while
// hydrostatic balance needs dp/dr = B there. That over-determines the discrete Neumann problem, and
// the same defect shows up as p_dyn drifting to +-14 bar at 300 sweeps instead of converging.
//
// THE SPLIT. p = p_hydro + p_dyn, with p_hydro defined by integration rather than by an elliptic
// solve:
//     p_hydro(r) = INTEGRAL from the base to r of the buoyancy,     dp_hydro/dr = buoyancy
// The radial force is then balanced exactly and by construction, and what remains of the buoyancy
// in the momentum equation is the HORIZONTAL gradient of p_hydro, which is the term that actually
// drives a circulation. On this grid that is about 0.012 against transport of order one, because
// the horizontal derivative carries a factor 1/r with r = 500 — small, baroclinic, and physically
// what a thin shell on a large planet should feel. p_dyn is left with the barotropic and
// non-hydrostatic remainder, which is precisely the part it CAN represent, because that part does
// not require a nonzero normal gradient at the walls.
//
// WHAT IS GIVEN UP. Exactly balancing the radial force makes the model hydrostatic in the vertical:
// buoyancy no longer accelerates u directly. At 3.5 km vertical against ~1000 km horizontal
// resolution that is the correct approximation and the one every large-scale GCM makes; a model
// that cannot solve for a non-hydrostatic pressure gains nothing by pretending to carry one.
//
// THE CONSTANT OF INTEGRATION is a free function of (theta, phi) and it is not cosmetic: it sets
// the barotropic part of the horizontal pressure gradient, i.e. a depth-independent horizontal
// force felt through the whole column. It is a modelling decision about where the atmosphere is
// anchored, and it has been taken: **the base is the reference isobaric surface**, so the integral
// runs upward from it. That is the natural anchor here — the base of this shell is the deep, dense
// boundary at 11 bar, where the gas is some fifty times denser than at the top and horizontal
// pressure contrasts are correspondingly harder to sustain, while the model top at 0.02 bar is an
// arbitrary cut through a continuing atmosphere with nothing to anchor it. ATJUP_HYDRO_REF=1
// integrates downward from the top instead, and exists so the choice can be measured rather than
// argued; it is not the intended configuration.
//
// Computed once per Runge-Kutta step, like buoy_ref_level, so it is the state at the start of the
// step and is held fixed through the four stages.
void cJupiterModel::computeHydrostaticPressure(){
    static const bool ref_top = [](){
        const char* e = getenv("ATJUP_HYDRO_REF"); return e && atoi(e) != 0; }();

    // The same factor rhs_u puts on the buoyancy, read the same way, so the two cannot drift apart.
    static const double nd_buoy_local = [](){
        const char* e = getenv("ATJUP_ND_BUOY");
        const char* n = getenv("ATJUP_NONDIM");
        const int on = e ? atoi(e) : (n ? atoi(n) : 0);
        return on != 0 ? 1.0 : 0.0; }();
    static const double buoy_scale_local = [](){
        const char* e = getenv("ATJUP_BUOY_SCALE"); return e ? atof(e) : 1.0; }();

    const double nd = (nd_buoy_local != 0.0) ? 1.0e5 * (L_atm * 1.0e3) / (u_0 * u_0) : 1.0;

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            // Only the fluid part of the column carries a hydrostatic integral; solid cells hold
            // bcSolidGround values and are left at zero so their gradients cannot leak into a
            // neighbouring fluid cell as a spurious force.
            int i0 = 0;
            while(i0 < im && SeaMount.x[i0][j][k] == 1.0){ p_hydro.x[i0][j][k] = 0.0; i0++; }
            if(i0 >= im) continue;
            int i1 = i0;
            while(i1 + 1 < im && SeaMount.x[i1 + 1][j][k] != 1.0) i1++;
            for(int i = i1 + 1; i < im; i++) p_hydro.x[i][j][k] = 0.0;

            auto buoy = [&](int i)->double{
                const double T = t.x[i][j][k];
                if(!(T > 0.0)) return 0.0;
                const double b = g * buoy_pressure(i, j, k) / (r_mix * R_mix * T * t_ref)
                               - buoy_ref_level[i];
                const double f = -nd * buoy_scale_local * buoyancy * b;
                return std::isfinite(f) ? f : 0.0;
            };

            if(!ref_top){
                p_hydro.x[i0][j][k] = 0.0;
                for(int i = i0 + 1; i <= i1; i++)
                    p_hydro.x[i][j][k] = p_hydro.x[i-1][j][k]
                                       + 0.5 * (buoy(i-1) + buoy(i)) * dr;
            } else {
                p_hydro.x[i1][j][k] = 0.0;
                for(int i = i1 - 1; i >= i0; i--)
                    p_hydro.x[i][j][k] = p_hydro.x[i+1][j][k]
                                       - 0.5 * (buoy(i+1) + buoy(i)) * dr;
            }
        }
    }
}

void cJupiterModel::RungeKuttaJup(){
    cout << endl << "      ATJUP: RungeKuttaJup" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    computeBuoyancyRefLevel();     // refresh the buoyancy base state for this RK4 step
    computeHydrostaticPressure();  // and the hydrostatic pressure built from it

    // ---- k* ceiling, as in ATOM's turbulent RK4 ----
    // The k production term is linear in k while its sink beta*.k.omega is linear too, so an
    // unbalanced production region can grow k exponentially; ATOM hit 1e98 m2/s2 at a single
    // polar cell before capping. The cap is a runaway guard, deliberately far above anything
    // physical: the closure's own Jovian equilibrium is k ~ 67 m2/s2 (see ParaView_Jup.cpp),
    // and 1000 m2/s2 is above the strongest hurricane-core TKE on record. Override with
    // ATJUP_TKE_MAX [m2/s2] if a storm study needs more headroom.
    static const double tke_max_phys = [](){
        const char* e = getenv("ATJUP_TKE_MAX"); return e ? atof(e) : 1000.0; }();
    const double tke_max_nd = tke_max_phys / (u_0 * u_0);
    constexpr double dis_min_nd = 1.0e-10;    // matches TurbulenceJup::dis_min

    // NaN-safe clamp. std::min/std::max compare with unguarded `<`; under -ffast-math
    // (-ffinite-math-only) the compiler may assume the operands are finite, so a NaN can slip
    // through a plain clamp. Detect non-finite values from the IEEE-754 exponent bits instead
    // and fall back to the lower bound. Same trick ATOM uses in RungeKutta_Atm_Turb.cpp.
    auto safe_clamp = [](double v, double lo, double hi) -> double {
        std::uint64_t bits;
        std::memcpy(&bits, &v, sizeof(bits));
        if((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) return lo;
        return (v < lo) ? lo : ((v > hi) ? hi : v);
    };

    // Precompute sin/cos tables — depend only on j
    // sinthe is clamped to a minimum to prevent 1/sin²θ blow-up near the poles.
    // The sequential k-loop in the RK creates an asymmetric phi Laplacian whose
    // error is amplified by inv_rm2sinthe2; clamping keeps the amplification < 1.
    // PressureSolverJup uses the same threshold. Raised 0.4 -> 0.55 (ATOM parity, metric
    // floor ~57°) to curb the polar 1/sin²θ amplification that seeded a long-run pole blow-up.
    const double sinthe_min = cJupiterModel::sinthe_min();   // see the note in cJupiterModel.h
    std::vector<double> sinthe_tbl(jm), costhe_tbl(jm);
    for(int j = 0; j < jm; j++){
        sinthe_tbl[j] = std::max(sinthe_min, std::abs(sin(the.z[j])));
        costhe_tbl[j] = cos(the.z[j]);
    }

    const double inv_2dr   = 1.0 / (2.0 * dr);
    const double inv_2dthe = 1.0 / (2.0 * dthe);
    const double inv_2dphi = 1.0 / (2.0 * dphi);
    const double inv_dr2   = 1.0 / (dr   * dr);
    const double inv_dthe2 = 1.0 / (dthe * dthe);
    const double inv_dphi2 = 1.0 / (dphi * dphi);

    static const int nh4sh_cap_on = [](){
        const char* e = getenv("ATJUP_NH4SH_CAP"); return e ? atoi(e) : 0; }();
    long nh4sh_cap_hits = 0;

    // ====================================================================================
    // Monotonicity limiter for the temperature transport (ATJUP_T_LIMITER, default 0 = off,
    // every existing run bit-identical).
    //
    // WHY. The advection of t is centred and unlimited, and centred differences are not
    // monotone: at a sharp front they under- and overshoot. Measured, dt=0.025, 200 sweeps,
    // level i=20, minimum of t over the level:
    //
    //     initial global minimum anywhere in the model   110.0000 K  (-163.1500 degC)
    //     without the obstacle, iter 55..100             -163.19 .. -163.30 degC, flat
    //     with the obstacle,     iter 60 / 65 / 85 / 90  -166.2 / -171.5 / -171.8 / -581842
    //
    // Without the cone the field transports the cold isothermal stratosphere down, lands on
    // its value and stays there — correct. With the cone it goes 8.4 K BELOW the coldest value
    // that exists anywhere in the initial data, and from there to NaN at iteration 89.
    // Advection cannot create a new extremum; a field that leaves the range of its own initial
    // data is a scheme error, and the sharpest front in the model is the staircase wall of the
    // cone. The 110 K itself is imposed exactly once, by the clamp in init_temperature
    // (InitValues_Jup.cpp:245) — there is no relaxation term anywhere in the time loop.
    //
    // What this does: after the RK4 update, clip t to the range spanned by the cell and its six
    // face neighbours IN THE OLD STATE (tn, untouched during the update). That is the clipping
    // step of an FCT scheme — it cannot create a new local extremum, and it leaves any update
    // that stays inside the local bounds exactly as it was. Solid neighbours are skipped: their
    // tn is an extrapolated value, not a state. It reports how often it bites, because a limiter
    // that acts silently would hide the front it is standing in for.
    //
    // The species already have their guard (FluxLimiterNH4SH, damp_wiggles on the mass fluxes);
    // t had none, and is also the only prognostic field with no bound of any kind.
    // ====================================================================================
    static const int t_limiter_on = [](){
        const char* e = getenv("ATJUP_T_LIMITER"); return e ? atoi(e) : 0; }();
    long t_clip_hits = 0;

    // ===== RK4, WITH THE FOUR STAGES SEPARATED =====
    //
    // Each stage is now two passes over the grid with an implicit barrier between them: one fills
    // rhs_* from the current stage input, one folds rhs_* into an accumulator and forms the next
    // stage input. Eight passes per step instead of one.
    //
    // WHY. All four stages used to run inside a single cell loop: RHSJup(i,j,k) then an immediate
    // overwrite of t, u, v, w and every species AT THAT CELL, four times over, before the loop
    // moved on. But RHSJup DIFFERENTIATES those same live fields at i+-1, j+-1, k+-1, so a cell
    // computing its stage 2 read neighbours that had already been advanced to their own stage 1,
    // or had not been touched yet, depending on where the loop had got to.
    //
    // Under OpenMP that is a data race: two runs of this binary at 24 threads differed in 5 of the
    // 7 output files. Serially it was not a race but it was still not RK4 — every stage
    // differentiated neighbours sitting at inconsistent stages, so the scheme was a pointwise
    // four-substep update wearing RK4's coefficients. One cause, one fix.
    //
    // Found first in ATSAT, whose loop is the same construction, by bisecting the model one
    // parallel region at a time; this file has the identical defect and now the identical repair.
    //
    // SOLID CELLS are skipped in both passes exactly as the single loop skipped them, so the
    // SeaMount interior is left untouched and no tendency is ever evaluated inside the body.
    for(int stage = 0; stage < 4; stage++){

        const double c_in = (stage == 0 || stage == 1) ? 0.5 * dt : (stage == 2 ? dt : 0.0);
        const double wgt  = (stage == 0 || stage == 3) ? 1.0 : 2.0;

        // ---- pass A: tendencies everywhere, from one consistent state ----
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < im-1; i++){
            for(int j = 3; j < jm-3; j++){
                CellGeometry geo;
                geo.rm      = rad.z[i];
                geo.rm2     = geo.rm * geo.rm;
                geo.exp_rm   = coord_stretching ? 1.0 / (geo.rm + 1.0) : 1.0;
                geo.exp_2_rm = geo.exp_rm * geo.exp_rm;
                geo.sinthe  = sinthe_tbl[j];
                geo.sinthe2 = geo.sinthe * geo.sinthe;
                geo.costhe  = costhe_tbl[j];
                if(costhe_abs() && j > 90) geo.costhe = -geo.costhe;
                geo.cotanthe            = geo.costhe / geo.sinthe;
                geo.inv_rm              = 1.0 / geo.rm;
                geo.inv_rm2             = 1.0 / geo.rm2;
                geo.inv_rmsinthe        = 1.0 / (geo.rm * geo.sinthe);
                geo.inv_rm2sinthe       = geo.inv_rm2 / geo.sinthe;
                geo.inv_rm2sinthe2      = geo.inv_rm2 / geo.sinthe2;
                geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;
                geo.inv_2dr   = inv_2dr;
                geo.inv_2dthe = inv_2dthe;
                geo.inv_2dphi = inv_2dphi;
                geo.inv_dr2   = inv_dr2;
                geo.inv_dthe2 = inv_dthe2;
                geo.inv_dphi2 = inv_dphi2;
                for(int k = 1; k < km-1; k++){
                    if(SeaMount.x[i][j][k] == 1.0) continue;
                    cJupiterModel::RHSJup(i, j, k, geo);
                }
            }
        }

        // ---- pass B: fold into the accumulator, then form the next stage input ----
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < im-1; i++){
            for(int j = 3; j < jm-3; j++){
                for(int k = 1; k < km-1; k++){
                    if(SeaMount.x[i][j][k] == 1.0) continue;
                    acc_t.x[i][j][k] = (stage == 0) ? wgt * rhs_t.x[i][j][k]
                                      : acc_t.x[i][j][k] + wgt * rhs_t.x[i][j][k];
                    acc_u.x[i][j][k] = (stage == 0) ? wgt * rhs_u.x[i][j][k]
                                      : acc_u.x[i][j][k] + wgt * rhs_u.x[i][j][k];
                    acc_v.x[i][j][k] = (stage == 0) ? wgt * rhs_v.x[i][j][k]
                                      : acc_v.x[i][j][k] + wgt * rhs_v.x[i][j][k];
                    acc_w.x[i][j][k] = (stage == 0) ? wgt * rhs_w.x[i][j][k]
                                      : acc_w.x[i][j][k] + wgt * rhs_w.x[i][j][k];
                    acc_h2o.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o.x[i][j][k]
                                      : acc_h2o.x[i][j][k] + wgt * rhs_h2o.x[i][j][k];
                    acc_h2o_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o_cloud.x[i][j][k]
                                      : acc_h2o_cloud.x[i][j][k] + wgt * rhs_h2o_cloud.x[i][j][k];
                    acc_h2o_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o_ice.x[i][j][k]
                                      : acc_h2o_ice.x[i][j][k] + wgt * rhs_h2o_ice.x[i][j][k];
                    acc_h2s.x[i][j][k] = (stage == 0) ? wgt * rhs_h2s.x[i][j][k]
                                      : acc_h2s.x[i][j][k] + wgt * rhs_h2s.x[i][j][k];
                    acc_nh3.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3.x[i][j][k]
                                      : acc_nh3.x[i][j][k] + wgt * rhs_nh3.x[i][j][k];
                    acc_nh3_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3_cloud.x[i][j][k]
                                      : acc_nh3_cloud.x[i][j][k] + wgt * rhs_nh3_cloud.x[i][j][k];
                    acc_nh3_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3_ice.x[i][j][k]
                                      : acc_nh3_ice.x[i][j][k] + wgt * rhs_nh3_ice.x[i][j][k];
                    acc_ch4.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4.x[i][j][k]
                                      : acc_ch4.x[i][j][k] + wgt * rhs_ch4.x[i][j][k];
                    acc_ch4_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4_cloud.x[i][j][k]
                                      : acc_ch4_cloud.x[i][j][k] + wgt * rhs_ch4_cloud.x[i][j][k];
                    acc_ch4_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4_ice.x[i][j][k]
                                      : acc_ch4_ice.x[i][j][k] + wgt * rhs_ch4_ice.x[i][j][k];
                    acc_nh4sh.x[i][j][k] = (stage == 0) ? wgt * rhs_nh4sh.x[i][j][k]
                                      : acc_nh4sh.x[i][j][k] + wgt * rhs_nh4sh.x[i][j][k];
                    acc_tke.x[i][j][k] = (stage == 0) ? wgt * rhs_tke.x[i][j][k]
                                      : acc_tke.x[i][j][k] + wgt * rhs_tke.x[i][j][k];
                    acc_dis.x[i][j][k] = (stage == 0) ? wgt * rhs_dis.x[i][j][k]
                                      : acc_dis.x[i][j][k] + wgt * rhs_dis.x[i][j][k];

                    if(stage < 3){
                        t.x[i][j][k] = tn.x[i][j][k] + c_in * rhs_t.x[i][j][k];
                        u.x[i][j][k] = un.x[i][j][k] + c_in * rhs_u.x[i][j][k];
                        v.x[i][j][k] = vn.x[i][j][k] + c_in * rhs_v.x[i][j][k];
                        w.x[i][j][k] = wn.x[i][j][k] + c_in * rhs_w.x[i][j][k];
                        h2o.x[i][j][k] = h2on.x[i][j][k] + c_in * rhs_h2o.x[i][j][k];
                        h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] + c_in * rhs_h2o_cloud.x[i][j][k];
                        h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] + c_in * rhs_h2o_ice.x[i][j][k];
                        h2s.x[i][j][k] = h2sn.x[i][j][k] + c_in * rhs_h2s.x[i][j][k];
                        nh3.x[i][j][k] = nh3n.x[i][j][k] + c_in * rhs_nh3.x[i][j][k];
                        nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] + c_in * rhs_nh3_cloud.x[i][j][k];
                        nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] + c_in * rhs_nh3_ice.x[i][j][k];
                        ch4.x[i][j][k] = ch4n.x[i][j][k] + c_in * rhs_ch4.x[i][j][k];
                        ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] + c_in * rhs_ch4_cloud.x[i][j][k];
                        ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] + c_in * rhs_ch4_ice.x[i][j][k];
                        nh4sh.x[i][j][k] = nh4shn.x[i][j][k] + c_in * rhs_nh4sh.x[i][j][k];
                        // k* is positive-definite and dis* strictly positive; clamp at every stage
                        // so a stiff intermediate can never feed a negative k or a zero omega back
                        // into the next RHS evaluation (nue* = k/omega would then be non-finite).
                        tke.x[i][j][k] = safe_clamp(tken.x[i][j][k]
                            + c_in * rhs_tke.x[i][j][k], 0.0, tke_max_nd);
                        dis.x[i][j][k] = std::max(dis_min_nd, disn.x[i][j][k]
                            + c_in * rhs_dis.x[i][j][k]);
                    }
                }
            }
        }
    }

    // ===== Final assembly: y_{n+1} = y_n + dt/6 (k1 + 2k2 + 2k3 + k4) =====
    {
        const double one_sixth = 1.0 / 6.0;
        #pragma omp parallel for collapse(2) schedule(static) \
                reduction(+:t_clip_hits, nh4sh_cap_hits)
        for(int i = 1; i < im-1; i++){
            for(int j = 3; j < jm-3; j++){
                for(int k = 1; k < km-1; k++){
                    if(SeaMount.x[i][j][k] == 1.0) continue;

                    const double tn_ijk = tn.x[i][j][k];
                    {
                        double t_new = tn_ijk + dt * acc_t.x[i][j][k] * one_sixth;
                        if(t_limiter_on){
                            double lo = tn_ijk, hi = tn_ijk;
                            auto take = [&](int ii, int jj, int kk){
                                const double v = tn.x[ii][jj][kk];
                                if(v < lo) lo = v;
                                if(v > hi) hi = v;
                            };
                            take(i-1,j,k); take(i+1,j,k);
                            take(i,j-1,k); take(i,j+1,k);
                            take(i,j,k-1); take(i,j,k+1);
                            if(t_new < lo){ t_new = lo; ++t_clip_hits; }
                            else if(t_new > hi){ t_new = hi; ++t_clip_hits; }
                        }
                        t.x[i][j][k] = t_new;
                    }
                    u.x[i][j][k] = un.x[i][j][k] + dt * acc_u.x[i][j][k] * one_sixth;
                    v.x[i][j][k] = vn.x[i][j][k] + dt * acc_v.x[i][j][k] * one_sixth;
                    w.x[i][j][k] = wn.x[i][j][k] + dt * acc_w.x[i][j][k] * one_sixth;
                    h2o.x[i][j][k] = std::max(0.0, h2on.x[i][j][k] + dt * acc_h2o.x[i][j][k] * one_sixth);
                    h2o_cloud.x[i][j][k] = std::max(0.0, h2o_cloudn.x[i][j][k] + dt * acc_h2o_cloud.x[i][j][k] * one_sixth);
                    h2o_ice.x[i][j][k] = std::max(0.0, h2o_icen.x[i][j][k] + dt * acc_h2o_ice.x[i][j][k] * one_sixth);
                    h2s.x[i][j][k] = std::max(0.0, h2sn.x[i][j][k] + dt * acc_h2s.x[i][j][k] * one_sixth);
                    nh3.x[i][j][k] = std::max(0.0, nh3n.x[i][j][k] + dt * acc_nh3.x[i][j][k] * one_sixth);
                    nh3_cloud.x[i][j][k] = std::max(0.0, nh3_cloudn.x[i][j][k] + dt * acc_nh3_cloud.x[i][j][k] * one_sixth);
                    nh3_ice.x[i][j][k] = std::max(0.0, nh3_icen.x[i][j][k] + dt * acc_nh3_ice.x[i][j][k] * one_sixth);
                    ch4.x[i][j][k] = std::max(0.0, ch4n.x[i][j][k] + dt * acc_ch4.x[i][j][k] * one_sixth);
                    ch4_cloud.x[i][j][k] = std::max(0.0, ch4_cloudn.x[i][j][k] + dt * acc_ch4_cloud.x[i][j][k] * one_sixth);
                    ch4_ice.x[i][j][k] = std::max(0.0, ch4_icen.x[i][j][k] + dt * acc_ch4_ice.x[i][j][k] * one_sixth);
                    // NH4SH ceiling (ATJUP_NH4SH_CAP, DEFAULT 0 = off). A dressing, not a cure; it
                    // belongs after the two rate-law repairs in ChemistryJup.h.
                    //
                    // IT IS OFF BECAUSE THE BOUND BELOW IS WRONG, and the measurement says so.
                    // The intent was r_max[j], the bound the INITIAL field is built with
                    // (InitValues_Jup.cpp:536). But `r_max` is a SCRATCH vector that every species
                    // initialiser rebuilds for itself — CH4, H2O, NH3 and finally H2S at
                    // InitValues_Jup.cpp:611, which is the one that survives into the time loop.
                    // By the time RungeKuttaJup reads it, r_max[j] carries the H2S bound, not an
                    // NH4SH bound. Switched on it therefore bit in 385 000 to 513 000 cells PER
                    // ITERATION and drove NH4SH to zero everywhere.
                    //
                    // A correct ceiling has to be built for NH4SH and kept (e.g. the maximum of
                    // the initial NH4SH field, stored once at init), which is a decision about the
                    // model rather than a line of code.
                    {
                        double nh4sh_new = nh4shn.x[i][j][k]
                                         + dt * acc_nh4sh.x[i][j][k] * one_sixth;
                        if(nh4sh_new < 0.0) nh4sh_new = 0.0;
                        if(nh4sh_cap_on && nh4sh_new > r_max[j]){
                            nh4sh_new = r_max[j];
                            ++nh4sh_cap_hits;
                        }
                        nh4sh.x[i][j][k] = nh4sh_new;
                    }
                    tke.x[i][j][k] = safe_clamp(tken.x[i][j][k]
                        + dt * acc_tke.x[i][j][k] * one_sixth, 0.0, tke_max_nd);
                    dis.x[i][j][k] = std::max(dis_min_nd, disn.x[i][j][k]
                        + dt * acc_dis.x[i][j][k] * one_sixth);
                }
            }
        }
    }

    if(nh4sh_cap_hits > 0)
        printf("      ATJUP: NH4SH cap bit in %ld cells this iteration (ATJUP_NH4SH_CAP=0 to lift it)\n",
               nh4sh_cap_hits);
    if(t_clip_hits > 0)
        printf("      ATJUP: t limiter clipped %ld cells this iteration (ATJUP_T_LIMITER=0 to lift it)\n",
               t_clip_hits);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for solveRungeKutta\n", elapsed.count() * 1e-9);

    cout << "      ATJUP: RungeKuttaJup ended" << endl;
    return;
}
/*
*
*/
