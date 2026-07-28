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

    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < im-1; i++){
        for(int j = 3; j < jm-3; j++){

            // Build geometry struct once per (i,j)
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

                double tn_ijk     = tn.x[i][j][k];
                double un_ijk     = un.x[i][j][k];
                double vn_ijk     = vn.x[i][j][k];
                double wn_ijk     = wn.x[i][j][k];
                double h2on_ijk   = h2on.x[i][j][k];
                double h2ocn_ijk  = h2o_cloudn.x[i][j][k];
                double h2oin_ijk  = h2o_icen.x[i][j][k];
                double h2sn_ijk   = h2sn.x[i][j][k];
                double nh3n_ijk   = nh3n.x[i][j][k];
                double nh3cn_ijk  = nh3_cloudn.x[i][j][k];
                double nh3in_ijk  = nh3_icen.x[i][j][k];
                double ch4n_ijk   = ch4n.x[i][j][k];
                double ch4cn_ijk  = ch4_cloudn.x[i][j][k];
                double ch4in_ijk  = ch4_icen.x[i][j][k];
                double nh4shn_ijk = nh4shn.x[i][j][k];
                double tken_ijk   = tken.x[i][j][k];
                double disn_ijk   = disn.x[i][j][k];

                // ----- RK stage 1 -----
                cJupiterModel::RHSJup(i, j, k, geo);
                double kt1     = rhs_t.x[i][j][k];
                double ku1     = rhs_u.x[i][j][k];
                double kv1     = rhs_v.x[i][j][k];
                double kw1     = rhs_w.x[i][j][k];
                double kc1     = rhs_h2o.x[i][j][k];
                double kcloud1 = rhs_h2o_cloud.x[i][j][k];
                double kice1   = rhs_h2o_ice.x[i][j][k];
                double kh2s1   = rhs_h2s.x[i][j][k];
                double knh31   = rhs_nh3.x[i][j][k];
                double knh3c1  = rhs_nh3_cloud.x[i][j][k];
                double knh3i1  = rhs_nh3_ice.x[i][j][k];
                double kch41   = rhs_ch4.x[i][j][k];
                double kch4c1  = rhs_ch4_cloud.x[i][j][k];
                double kch4i1  = rhs_ch4_ice.x[i][j][k];
                double knh4sh1 = rhs_nh4sh.x[i][j][k];
                double ktke1   = rhs_tke.x[i][j][k];
                double kdis1   = rhs_dis.x[i][j][k];

                t.x[i][j][k]         = tn_ijk    + kt1     * 0.5 * dt;
                u.x[i][j][k]         = un_ijk    + ku1     * 0.5 * dt;
                v.x[i][j][k]         = vn_ijk    + kv1     * 0.5 * dt;
                w.x[i][j][k]         = wn_ijk    + kw1     * 0.5 * dt;
                h2o.x[i][j][k]       = h2on_ijk  + kc1     * 0.5 * dt;
                h2o_cloud.x[i][j][k] = h2ocn_ijk + kcloud1 * 0.5 * dt;
                h2o_ice.x[i][j][k]   = h2oin_ijk + kice1   * 0.5 * dt;
                h2s.x[i][j][k]       = h2sn_ijk  + kh2s1   * 0.5 * dt;
                nh3.x[i][j][k]       = nh3n_ijk  + knh31   * 0.5 * dt;
                nh3_cloud.x[i][j][k] = nh3cn_ijk + knh3c1  * 0.5 * dt;
                nh3_ice.x[i][j][k]   = nh3in_ijk + knh3i1  * 0.5 * dt;
                ch4.x[i][j][k]       = ch4n_ijk  + kch41   * 0.5 * dt;
                ch4_cloud.x[i][j][k] = ch4cn_ijk + kch4c1  * 0.5 * dt;
                ch4_ice.x[i][j][k]   = ch4in_ijk + kch4i1  * 0.5 * dt;
                nh4sh.x[i][j][k]     = nh4shn_ijk + knh4sh1 * 0.5 * dt;
                // k* is positive-definite and dis* strictly positive; clamp at every stage so a
                // stiff intermediate can never feed a negative k or a zero omega back into the
                // next RHS evaluation (nue* = k/omega would then be non-finite).
                tke.x[i][j][k]       = safe_clamp(tken_ijk + ktke1 * 0.5 * dt, 0.0, tke_max_nd);
                dis.x[i][j][k]       = std::max(dis_min_nd, disn_ijk + kdis1 * 0.5 * dt);

                // ----- RK stage 2 -----
                cJupiterModel::RHSJup(i, j, k, geo);
                double kt2     = rhs_t.x[i][j][k];
                double ku2     = rhs_u.x[i][j][k];
                double kv2     = rhs_v.x[i][j][k];
                double kw2     = rhs_w.x[i][j][k];
                double kc2     = rhs_h2o.x[i][j][k];
                double kcloud2 = rhs_h2o_cloud.x[i][j][k];
                double kice2   = rhs_h2o_ice.x[i][j][k];
                double kh2s2   = rhs_h2s.x[i][j][k];
                double knh32   = rhs_nh3.x[i][j][k];
                double knh3c2  = rhs_nh3_cloud.x[i][j][k];
                double knh3i2  = rhs_nh3_ice.x[i][j][k];
                double kch42   = rhs_ch4.x[i][j][k];
                double kch4c2  = rhs_ch4_cloud.x[i][j][k];
                double kch4i2  = rhs_ch4_ice.x[i][j][k];
                double knh4sh2 = rhs_nh4sh.x[i][j][k];
                double ktke2   = rhs_tke.x[i][j][k];
                double kdis2   = rhs_dis.x[i][j][k];

                t.x[i][j][k]         = tn_ijk    + kt2     * 0.5 * dt;
                u.x[i][j][k]         = un_ijk    + ku2     * 0.5 * dt;
                v.x[i][j][k]         = vn_ijk    + kv2     * 0.5 * dt;
                w.x[i][j][k]         = wn_ijk    + kw2     * 0.5 * dt;
                h2o.x[i][j][k]       = h2on_ijk  + kc2     * 0.5 * dt;
                h2o_cloud.x[i][j][k] = h2ocn_ijk + kcloud2 * 0.5 * dt;
                h2o_ice.x[i][j][k]   = h2oin_ijk + kice2   * 0.5 * dt;
                h2s.x[i][j][k]       = h2sn_ijk  + kh2s2   * 0.5 * dt;
                nh3.x[i][j][k]       = nh3n_ijk  + knh32   * 0.5 * dt;
                nh3_cloud.x[i][j][k] = nh3cn_ijk + knh3c2  * 0.5 * dt;
                nh3_ice.x[i][j][k]   = nh3in_ijk + knh3i2  * 0.5 * dt;
                ch4.x[i][j][k]       = ch4n_ijk  + kch42   * 0.5 * dt;
                ch4_cloud.x[i][j][k] = ch4cn_ijk + kch4c2  * 0.5 * dt;
                ch4_ice.x[i][j][k]   = ch4in_ijk + kch4i2  * 0.5 * dt;
                nh4sh.x[i][j][k]     = nh4shn_ijk + knh4sh2 * 0.5 * dt;
                tke.x[i][j][k]       = safe_clamp(tken_ijk + ktke2 * 0.5 * dt, 0.0, tke_max_nd);
                dis.x[i][j][k]       = std::max(dis_min_nd, disn_ijk + kdis2 * 0.5 * dt);

                // ----- RK stage 3 -----
                cJupiterModel::RHSJup(i, j, k, geo);
                double kt3     = rhs_t.x[i][j][k];
                double ku3     = rhs_u.x[i][j][k];
                double kv3     = rhs_v.x[i][j][k];
                double kw3     = rhs_w.x[i][j][k];
                double kc3     = rhs_h2o.x[i][j][k];
                double kcloud3 = rhs_h2o_cloud.x[i][j][k];
                double kice3   = rhs_h2o_ice.x[i][j][k];
                double kh2s3   = rhs_h2s.x[i][j][k];
                double knh33   = rhs_nh3.x[i][j][k];
                double knh3c3  = rhs_nh3_cloud.x[i][j][k];
                double knh3i3  = rhs_nh3_ice.x[i][j][k];
                double kch43   = rhs_ch4.x[i][j][k];
                double kch4c3  = rhs_ch4_cloud.x[i][j][k];
                double kch4i3  = rhs_ch4_ice.x[i][j][k];
                double knh4sh3 = rhs_nh4sh.x[i][j][k];
                double ktke3   = rhs_tke.x[i][j][k];
                double kdis3   = rhs_dis.x[i][j][k];

                t.x[i][j][k]         = tn_ijk    + kt3     * dt;
                u.x[i][j][k]         = un_ijk    + ku3     * dt;
                v.x[i][j][k]         = vn_ijk    + kv3     * dt;
                w.x[i][j][k]         = wn_ijk    + kw3     * dt;
                h2o.x[i][j][k]       = h2on_ijk  + kc3     * dt;
                h2o_cloud.x[i][j][k] = h2ocn_ijk + kcloud3 * dt;
                h2o_ice.x[i][j][k]   = h2oin_ijk + kice3   * dt;
                h2s.x[i][j][k]       = h2sn_ijk  + kh2s3   * dt;
                nh3.x[i][j][k]       = nh3n_ijk  + knh33   * dt;
                nh3_cloud.x[i][j][k] = nh3cn_ijk + knh3c3  * dt;
                nh3_ice.x[i][j][k]   = nh3in_ijk + knh3i3  * dt;
                ch4.x[i][j][k]       = ch4n_ijk  + kch43   * dt;
                ch4_cloud.x[i][j][k] = ch4cn_ijk + kch4c3  * dt;
                ch4_ice.x[i][j][k]   = ch4in_ijk + kch4i3  * dt;
                nh4sh.x[i][j][k]     = nh4shn_ijk + knh4sh3 * dt;
                tke.x[i][j][k]       = safe_clamp(tken_ijk + ktke3 * dt, 0.0, tke_max_nd);
                dis.x[i][j][k]       = std::max(dis_min_nd, disn_ijk + kdis3 * dt);

                // ----- RK stage 4 -----
                cJupiterModel::RHSJup(i, j, k, geo);
                double kt4     = rhs_t.x[i][j][k];
                double ku4     = rhs_u.x[i][j][k];
                double kv4     = rhs_v.x[i][j][k];
                double kw4     = rhs_w.x[i][j][k];
                double kc4     = rhs_h2o.x[i][j][k];
                double kcloud4 = rhs_h2o_cloud.x[i][j][k];
                double kice4   = rhs_h2o_ice.x[i][j][k];
                double kh2s4   = rhs_h2s.x[i][j][k];
                double knh34   = rhs_nh3.x[i][j][k];
                double knh3c4  = rhs_nh3_cloud.x[i][j][k];
                double knh3i4  = rhs_nh3_ice.x[i][j][k];
                double kch44   = rhs_ch4.x[i][j][k];
                double kch4c4  = rhs_ch4_cloud.x[i][j][k];
                double kch4i4  = rhs_ch4_ice.x[i][j][k];
                double knh4sh4 = rhs_nh4sh.x[i][j][k];
                double ktke4   = rhs_tke.x[i][j][k];
                double kdis4   = rhs_dis.x[i][j][k];

                // ----- Final RK4 update -----
                const double one_sixth = 1.0 / 6.0;
                t.x[i][j][k]         = tn_ijk    + dt * (kt1     + 2.0*kt2     + 2.0*kt3     + kt4    ) * one_sixth;
                u.x[i][j][k]         = un_ijk    + dt * (ku1     + 2.0*ku2     + 2.0*ku3     + ku4    ) * one_sixth;
                v.x[i][j][k]         = vn_ijk    + dt * (kv1     + 2.0*kv2     + 2.0*kv3     + kv4    ) * one_sixth;
                w.x[i][j][k]         = wn_ijk    + dt * (kw1     + 2.0*kw2     + 2.0*kw3     + kw4    ) * one_sixth;
                h2o.x[i][j][k]       = std::max(0.0, h2on_ijk  + dt * (kc1     + 2.0*kc2     + 2.0*kc3     + kc4    ) * one_sixth);
                h2o_cloud.x[i][j][k] = std::max(0.0, h2ocn_ijk + dt * (kcloud1 + 2.0*kcloud2 + 2.0*kcloud3 + kcloud4) * one_sixth);
                h2o_ice.x[i][j][k]   = std::max(0.0, h2oin_ijk + dt * (kice1   + 2.0*kice2   + 2.0*kice3   + kice4  ) * one_sixth);
                h2s.x[i][j][k]       = std::max(0.0, h2sn_ijk  + dt * (kh2s1   + 2.0*kh2s2   + 2.0*kh2s3   + kh2s4  ) * one_sixth);
                nh3.x[i][j][k]       = std::max(0.0, nh3n_ijk  + dt * (knh31   + 2.0*knh32   + 2.0*knh33   + knh34  ) * one_sixth);
                nh3_cloud.x[i][j][k] = std::max(0.0, nh3cn_ijk + dt * (knh3c1  + 2.0*knh3c2  + 2.0*knh3c3  + knh3c4 ) * one_sixth);
                nh3_ice.x[i][j][k]   = std::max(0.0, nh3in_ijk + dt * (knh3i1  + 2.0*knh3i2  + 2.0*knh3i3  + knh3i4 ) * one_sixth);
                ch4.x[i][j][k]       = std::max(0.0, ch4n_ijk  + dt * (kch41   + 2.0*kch42   + 2.0*kch43   + kch44  ) * one_sixth);
                ch4_cloud.x[i][j][k] = std::max(0.0, ch4cn_ijk + dt * (kch4c1  + 2.0*kch4c2  + 2.0*kch4c3  + kch4c4 ) * one_sixth);
                ch4_ice.x[i][j][k]   = std::max(0.0, ch4in_ijk + dt * (kch4i1  + 2.0*kch4i2  + 2.0*kch4i3  + kch4i4 ) * one_sixth);
                nh4sh.x[i][j][k]     = std::max(0.0, nh4shn_ijk + dt * (knh4sh1 + 2.0*knh4sh2 + 2.0*knh4sh3 + knh4sh4) * one_sixth);
                tke.x[i][j][k]       = safe_clamp(tken_ijk + dt * (ktke1 + 2.0*ktke2 + 2.0*ktke3 + ktke4) * one_sixth,
                                                  0.0, tke_max_nd);
                dis.x[i][j][k]       = std::max(dis_min_nd,
                                                disn_ijk + dt * (kdis1 + 2.0*kdis2 + 2.0*kdis3 + kdis4) * one_sixth);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for solveRungeKutta\n", elapsed.count() * 1e-9);

    cout << "      ATJUP: RungeKuttaJup ended" << endl;
    return;
}
/*
*
*/
