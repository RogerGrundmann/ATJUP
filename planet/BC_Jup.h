/*
 * Atmosphere General Circulation Modell (ATJUP)
 * Standalone boundary-condition class for the Jupiter model.
 * Declared as friend of cJupiterModel so it may access all private members
 * through the stored reference.
 *
 * This is a header-only file: it contains the BC_Jup class, all inline method
 * bodies, and the inline cJupiterModel delegation wrappers that forward each
 * cJupiterModel::BC_xxx() call to the corresponding BC_Jup method.
 * BC_Jup.cpp is a minimal stub that provides a translation unit.
*/

#pragma once

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <iostream>

class cJupiterModel;
class Array;

using namespace std;

class BC_Jup {
public:

    explicit BC_Jup(cJupiterModel& model) : m(model) {}

    void bcRadius();
    void bcTheta();
    void bcPhi();
    void bcSeaMount();
    void bcSolidGround();
    void bcVelSurfSur();
    void bcScalarSurfSur();
    void initTropopauseLayers();

private:
    cJupiterModel& m;
};


// -----------------------------------------------------------------------
// Inline implementation  (header-only, like ATOM's BC_Atm.h)
// -----------------------------------------------------------------------
#include "cJupiterModel.h"
#include <cstdlib>   // getenv/atoi for the boundary-hardening knobs

// ---------------------------------------------------------------------------
// Boundary-hardening knobs, ported from ATOM_Precipitation's BC_Atm.h after a
// side-by-side comparison. Set any knob to 0 to A/B against the old behaviour.
//
// Defaults reflect how well each is actually supported in ATJUP:
//   rigid_lid  1  — required for well-posedness (closed shell + all-Neumann Poisson),
//                   though measured effect at 100 iters was negligible.
//   top_taper  1  — MEASURED: moves the w maximum off the lid (133 -> 126 km) and
//                   restores the IC's intent. The one clear win of the four.
//   pole_copy  1  — numerically strictly better at the sin(theta) -> 0 singularity;
//                   no ATJUP symptom observed yet, kept as insurance.
//   seam_damp  0  — no ATJUP evidence of a seam mode; see its note below.
//
// A NOTE ON EVIDENCE, so nobody repeats the mistake: the clustering of printMinMax
// extrema at "0/1 deg E" and "90 deg N" is NOT evidence of boundary trouble. It is
// a tie-breaking artifact — searchMinMax_3D uses a strict > and scans k before i,
// so for any zonally near-uniform field (which the IC is, by construction) the max
// is reported at the smallest k that attains it, i.e. k = 0 or 1, and at j = 0.
// Verified: enabling the pole and seam treatments below changed that clustering not
// at all (28 -> 29 seam, 17 -> 17 polar out of 46). Judge these knobs with direct
// field diagnostics, not with extremum locations.
// ---------------------------------------------------------------------------
namespace BCJupKnobs {
    inline int env_int(const char* name, int dflt){
        const char* e = getenv(name);
        return e ? atoi(e) : dflt;
    }
    // (1) u = 0 at i=0 and i=im-1 (no flow through either radial wall).
    inline int rigid_lid()  { static const int v = env_int("ATJUP_BC_RIGID_LID", 1); return v; }
    // (2) taper v,w to zero over the top three layers.
    inline int top_taper()  { static const int v = env_int("ATJUP_BC_TOP_TAPER", 1); return v; }
    // (3) plain copy instead of (4/3,-1/3) extrapolation at the poles.
    inline int pole_copy()  { static const int v = env_int("ATJUP_BC_POLE_COPY", 1); return v; }
    // (7) Radial boundary: plain copy f[0]=f[1], f[im-1]=f[im-2] instead of the (4/3,-1/3)
    // extrapolation. DIAGNOSTIC KNOB, default 0 (off). The two agree for a flat profile, but on
    // a profile that varies across the wall the 2-point form overshoots by 1/3 of the interior
    // increment, and the deep-level w growth that ends the long runs sits exactly at i = 0..2.
    // Running with 1 turns the overshoot off without touching anything else, which separates a
    // boundary artifact from a genuine missing momentum sink. It is deliberately NOT the default:
    // the copy is only first-order accurate and would degrade every transported field.
    inline int radius_copy(){ static const int v = env_int("ATJUP_BC_RADIUS_COPY", 0); return v; }

    // (8) Obstacle surface: plain copy f[s] = f[a] instead of the (4/3,-1/3) extrapolation for
    // the scalars in the SeaMount surface cells. DIAGNOSTIC KNOB, default 0 (off).
    //
    // WHY IT EXISTS. The fluid cells flanking the cone cool without stopping. Measured on level
    // i=20, both flanks, pairwise symmetric about the cone longitude k_0=180:
    //
    //   dt=0.001, 450 iter (10.5 min Jupiter time):  (95,166) -12.6 K   (105,159) -22.0 K
    //                                                (95,194) -12.5 K   (105,201) -22.0 K
    //   dt=0.025, 100 iter (52 min):   (94,166) -83.2 -> -153.6 degC, then NaN at iteration 89
    //
    // The rate per unit PHYSICAL time agrees between the two timesteps (4.08 K in 630 s against
    // 5.1 K in 700 s at the same cell), so it is not a timestep artifact — the large step only
    // reaches the end sooner. ATJUP_NANCHECK names t in a solid cell (20,94,170) as the first
    // non-finite value, before any velocity.
    //
    // The suspected loop runs across the wall: a fluid cell beside the obstacle still forms its
    // advection with CENTRED differences, so it reads the solid neighbour, whose value is this
    // extrapolation. With (4/3,-1/3) that value OVERSHOOTS past the near fluid cell whenever the
    // profile falls towards the wall, the fluid cell advects the too-cold ghost in, cools, and
    // the next extrapolation overshoots further. The comment at the lambda puts the per-step
    // amplification at 5/3 and expects fluid diffusion to bound it; on these timescales it does
    // not. A plain copy cannot overshoot, so this knob separates "the extrapolation amplifies it"
    // from "the advection scheme itself does" WITHOUT changing the advection.
    //
    // It is deliberately NOT the repair and not the default: the copy is only first-order and
    // the physical fix is one-sided differencing at the obstacle face, so that nothing is
    // differenced through the boundary at all.
    inline int mount_copy() { static const int v = env_int("ATJUP_BC_MOUNT_COPY", 0); return v; }

    // (4) number of 1-2-1 Shapiro passes across the phi seam. DEFAULT 0 (off): unlike the
    // others this one is NOT supported by ATJUP evidence. A direct zonal-roughness measurement
    // at iteration 100 found the seam SMOOTHER than the interior wherever real dynamics exist
    // (mean |d2w/dphi2| ratio seam:interior = 0.01-0.04 at 0, -22, -45 deg); it only stands out
    // at 45 deg N, where the field is nearly zonally uniform and the absolute amplitude is
    // ~1e-4 m/s. So there is no live seam instability to damp here, and leaving a filter that
    // perturbs u,v,w every iteration switched on would only confound later experiments.
    // Set to 2 if a long run ever shows zonal energy accumulating at k = 0/1/km-2.
    inline int seam_damp()  { static const int v = env_int("ATJUP_BC_SEAM_DAMP", 0); return v; }

    // (6) Lid temperature pin. DEFAULT 0 (off) — and the reason is worth recording, because the
    // motivation for porting this from ATOM turned out not to apply to ATJUP.
    //
    // ATOM pins t at i=im-1 because its cubic lid extrapolation projected interior curvature
    // onto the lid and its stratospheric top drifted upward within ~20 iterations, corrupting an
    // otherwise steady initial state. MEASURED in ATJUP instead (j=87, k=180): the lid is 56.34 K
    // in the initial condition and 57.53 K after 100 iterations — a drift of +1.2 K, negligible,
    // and in the direction we WANT. ATJUP's lid is not drifting.
    //
    // So the model's anomalously cold top (~57 K against a real 110-140 K) is NOT a boundary
    // artifact: it is the initial condition. init_temperature (InitValues_Jup.cpp) applies an
    // unbounded linear lapse rate, t = -gam*height/t_ref + t[0], with gam = 2.0 K/km over the
    // full L_atm = 140 km, i.e. T(lid) = T(0) - 280 K. There is no tropopause floor and no
    // stratospheric inversion, even though initTropopauseLayers() below computes a tropopause at
    // 125 km (equator) / 115 km (pole) that VelocityInitializerJup already honours when it ramps
    // v,w. The temperature IC simply ignores it.
    //
    // Pinning to the IC snapshot would therefore FREEZE the top at 56 K and stop the radiation
    // coupling from ever warming it — the opposite of what is wanted. Two useful modes are
    // nevertheless provided:
    //   ATJUP_BC_T_LID_PIN=1                  hold the lid at its IC value (ATOM parity; use to
    //                                         isolate lid drift, not to fix the cold top).
    //   ATJUP_BC_T_LID_PIN=1 ATJUP_BC_T_LID_K=115
    //                                         hold the lid at a PRESCRIBED physical temperature
    //                                         in kelvin, which is the version that actually
    //                                         serves the warm-top goal.
    inline int    t_lid_pin() { static const int    v = env_int("ATJUP_BC_T_LID_PIN", 0); return v; }
    inline double t_lid_K()   { static const double v = [](){ const char* e = getenv("ATJUP_BC_T_LID_K"); return e ? atof(e) : 0.0; }(); return v; }
}


inline void BC_Jup::bcRadius()
{
    const int im = m.im, jm = m.jm, km = m.km;
    const double c43 = m.c43, c13 = m.c13;

    Array* fields[] = {
        &m.t, &m.u, &m.v, &m.w,
        &m.ch4, &m.ch4_cloud, &m.ch4_ice,
        &m.h2o, &m.h2o_cloud, &m.h2o_ice,
        &m.h2s,
        &m.nh3, &m.nh3_cloud, &m.nh3_ice,
        &m.nh4sh,
        &m.j_h2s,  &m.j_nh3,  &m.j_nh4sh,
        &m.jT_h2s, &m.jT_nh3, &m.jT_nh4sh,
        &m.w_h2s,  &m.w_nh3,  &m.w_nh4sh,
        &m.massflux_h2s,  &m.massflux_nh3,  &m.massflux_nh4sh,
        &m.difflux_h2s,   &m.difflux_nh3,   &m.difflux_nh4sh,
        &m.thermalmassflux,
        &m.CoriolisForce, &m.CentrifugalForce,
        &m.BuoyancyForce, &m.PresGradForce,
        &m.Q_Latent, &m.Q_Sensible,
        // k* and dis* are prognostic (RHS_Jup_Turb.cpp), so their domain boundaries must be
        // extrapolated like every other transported field — ATOM carries them in the same
        // lists. Without this they stay frozen at the array initial values (k*=0, dis*=1e-10)
        // while the interior runs at dis* ~ 16, and the resulting permanent Laplacian at
        // i=0/im-1, the poles and the phi seam drains dis* to its floor within ~20 iterations.
        // nue* = k*/omega* then saturates its ceiling and, with ATJUP_TURB_COUPLING on, the
        // eddy viscosity destroys the momentum field.
        &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
    };
    const int nf = (int)(sizeof(fields) / sizeof(fields[0]));

    // 2-point Neumann extrapolation: f[s] = (4/3)f[a] - (1/3)f[b].
    // The 3-point cubic (3f[a]-3f[b]+f[c]) amplifies alternating errors by 7x
    // per call and blows up near the SeaMount contour (same reason bcSolidGround
    // was switched to the 2-point formula; bcRadius had the same latent bug).
    const int iml = im - 1;
    const bool do_lid   = BCJupKnobs::rigid_lid()  != 0;
    const bool do_taper = BCJupKnobs::top_taper()  != 0;
    const bool do_copy  = BCJupKnobs::radius_copy() != 0;   // (7) diagnostic: drop the overshoot

    // --- (6) Lid temperature pin (opt-in; see the knob note above for why it is off by
    // default and why ATJUP's cold top is an IC problem, not a boundary one). ---
    // The snapshot is taken on the FIRST call, before the extrapolation below overwrites
    // t.x[iml]: the first bcRadius() runs after all initialisation, so this captures the
    // initial condition. Done serially, outside the parallel region.
    const bool   do_t_pin = BCJupKnobs::t_lid_pin() != 0;
    const double t_lid_K  = BCJupKnobs::t_lid_K();
    if(do_t_pin && (int)m.t_top_init.size() != jm){
        m.t_top_init.assign(jm, std::vector<double>(km, 0.0));
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++)
                m.t_top_init[j][k] = (t_lid_K > 0.0) ? (t_lid_K / m.t_ref)
                                                     : m.t.x[iml][j][k];
    }
    const bool pin_t_top = do_t_pin && ((int)m.t_top_init.size() == jm);

    #pragma omp parallel for schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                if(do_copy){
                    F.x[0][j][k]    = F.x[1][j][k];
                    F.x[im-1][j][k] = F.x[im-2][j][k];
                }else{
                    F.x[0][j][k]    = c43*F.x[1][j][k]    - c13*F.x[2][j][k];
                    F.x[im-1][j][k] = c43*F.x[im-2][j][k] - c13*F.x[im-3][j][k];
                }
            }

            // --- (1) Rigid walls on the RADIAL velocity u at both radial boundaries. ---
            // u is the wall-normal component there, and the modelled shell is closed: no mass
            // crosses the lid, and none crosses the deep bottom either (RadiationJup injects the
            // interior heat flux F_int as ENERGY at i=0, not as mass). The Neumann extrapolation
            // just written is correct only for the TANGENTIAL v,w — applied to u it overshoots an
            // increasing radial profile and feeds back through the i=im-2 d/dr stencil, ratcheting
            // the vertical velocity up every step (ATOM measured u: 0 -> 17 m/s by iter 300 before
            // adding this; ATJUP currently reaches |u| ~ 33-38 m/s at 73-94 km, which is enormous
            // for a radial velocity in a 140 km shell).
            //
            // It also makes the pressure problem well posed: PressureSolverJup is all-Neumann, so
            // without u = 0 at both ends the column-mean vertical velocity is an undetermined,
            // freely drifting constant. Only u is pinned — v,w keep their Neumann (mirror) values,
            // which is the free-slip tangential condition.
            //
            // un is NOT written here: restoreVar(1.0) copies u -> un after all BCs run
            // (cJupiterModel.cpp), so the next RK4 step already starts from the wall value.
            if(do_lid){
                m.u.x[0][j][k]   = 0.0;
                m.u.x[iml][j][k] = 0.0;
            }

            // (6) Override the lid temperature extrapolation with the pinned value.
            if(pin_t_top) m.t.x[iml][j][k] = m.t_top_init[j][k];

            // --- (2) Taper the HORIZONTAL velocities to a quiet grid ceiling. ---
            // VelocityInitializerJup ramps v,w linearly to zero between the tropopause and the
            // model top, but the Neumann extrapolation above copies the interior value straight
            // to i=im-1 and so drags the zonal jet back up to the lid, undoing the IC ("stretched
            // up" jets). Symptom: ATJUP's w maximum sits at 133 km, i.e. i=38 of 40, right at the
            // lid where the IC intended zero. Ramping the top three layers by 2/3, 1/3, 0 restores
            // a quiet lid without the one-cell shear shock a hard zero would create.
            //
            // No runaway compounding: RK4 integrates i=1..im-2, so v,w at iml-1 and iml-2 are
            // recomputed from tendencies every iteration and the factor is re-applied to a fresh
            // value rather than to an already-tapered one. i=iml is outside the RK4 range, so
            // setting it to zero is a clean Dirichlet condition.
            if(do_taper){
                m.v.x[iml][j][k]    = 0.0;
                m.w.x[iml][j][k]    = 0.0;
                m.v.x[iml-1][j][k] *= (1.0/3.0);
                m.w.x[iml-1][j][k] *= (1.0/3.0);
                m.v.x[iml-2][j][k] *= (2.0/3.0);
                m.w.x[iml-2][j][k] *= (2.0/3.0);
            }
        }
    }
}


inline void BC_Jup::bcTheta()
{
    const int im = m.im, jm = m.jm, km = m.km;
    const double c43 = m.c43, c13 = m.c13;

    // All fields except v and w receive 2-point Neumann extrapolation at poles.
    Array* extrap_fields[] = {
        &m.t, &m.u,
        &m.ch4, &m.ch4_cloud, &m.ch4_ice,
        &m.h2o, &m.h2o_cloud, &m.h2o_ice,
        &m.h2s,
        &m.nh3, &m.nh3_cloud, &m.nh3_ice,
        &m.nh4sh,
        &m.j_h2s,  &m.j_nh3,  &m.j_nh4sh,
        &m.jT_h2s, &m.jT_nh3, &m.jT_nh4sh,
        &m.w_h2s,  &m.w_nh3,  &m.w_nh4sh,
        &m.thermalmassflux,
        &m.CoriolisForce, &m.CentrifugalForce,
        &m.BuoyancyForce, &m.PresGradForce,
        &m.Q_Latent, &m.Q_Sensible,
        // prognostic turbulence fields — see the note in bcRadius()
        &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
    };
    const int nf = (int)(sizeof(extrap_fields) / sizeof(extrap_fields[0]));

    // Flux fields that are singular-prone near poles: zero at pole boundary.
    Array* zero_at_poles[] = {
        &m.massflux_h2s, &m.massflux_nh3, &m.massflux_nh4sh,
        &m.fluxlim_nh4sh,
        &m.difflux_h2s,  &m.difflux_nh3,  &m.difflux_nh4sh,
    };
    const int nz = (int)(sizeof(zero_at_poles) / sizeof(zero_at_poles[0]));
    const bool pole_copy = BCJupKnobs::pole_copy() != 0;

    #pragma omp parallel for schedule(static)
//    for(int k = 1; k < km-1; k++){
//        for(int i = 1; i < im-1; i++){
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            m.v.x[i][0][k]    = 0.0;
            m.v.x[i][jm-1][k] = 0.0;
            m.w.x[i][0][k]    = 0.0;
            m.w.x[i][jm-1][k] = 0.0;

            // --- (3) Pole boundary: plain copy, not (4/3,-1/3) extrapolation. ---
            // At the spherical singularity sin(theta) -> 0 any extrapolation amplifies grid
            // noise: the (4/3,-1/3) form by 4/3 per call, the 3-point cubic by ~7. With a
            // strong bulk flow that is washed out, but with a weak or spinning-up velocity
            // field the amplified noise dominates — ATOM reached NaN at the pole within ~150
            // iterations before switching to a plain copy. A copy has amplification exactly
            // 1.0 and is the axisymmetric-pole assumption to first order, which is what the
            // pole physically is. Note this is justified on numerical grounds only — no ATJUP
            // diagnostic currently shows polar noise growth at 100 iterations, so treat it as
            // insurance for long / weakly-forced runs rather than a fix for an observed problem.
            if(pole_copy){
                for(int f = 0; f < nf; f++){
                    Array& F = *extrap_fields[f];
                    F.x[i][0][k]    = F.x[i][1][k];
                    F.x[i][jm-1][k] = F.x[i][jm-2][k];
                }
            } else {
                for(int f = 0; f < nf; f++){
                    Array& F = *extrap_fields[f];
                    F.x[i][0][k]    = c43*F.x[i][1][k]    - c13*F.x[i][2][k];
                    F.x[i][jm-1][k] = c43*F.x[i][jm-2][k] - c13*F.x[i][jm-3][k];
                }
            }

            for(int f = 0; f < nz; f++){
                Array& F = *zero_at_poles[f];
                F.x[i][0][k]    = 0.0;
                F.x[i][jm-1][k] = 0.0;
            }
        }
    }
}


inline void BC_Jup::bcPhi()
{
    const int im = m.im, jm = m.jm, km = m.km;
    const double c43 = m.c43, c13 = m.c13;

    Array* fields[] = {
        &m.t, &m.u, &m.v, &m.w,
        &m.ch4, &m.ch4_cloud, &m.ch4_ice,
        &m.h2o, &m.h2o_cloud, &m.h2o_ice,
        &m.h2s,
        &m.nh3, &m.nh3_cloud, &m.nh3_ice,
        &m.nh4sh,
        &m.j_h2s,  &m.j_nh3,
        &m.jT_h2s, &m.jT_nh3,
        &m.w_h2s,  &m.w_nh3,  &m.w_nh4sh,
        &m.massflux_h2s,  &m.massflux_nh3,  &m.massflux_nh4sh,
        &m.fluxlim_nh4sh,
        &m.difflux_h2s,   &m.difflux_nh3,   &m.difflux_nh4sh,
        &m.thermalmassflux,
        &m.CoriolisForce, &m.CentrifugalForce,
        &m.BuoyancyForce, &m.PresGradForce,
        &m.Q_Latent, &m.Q_Sensible,
        // prognostic turbulence fields — see the note in bcRadius()
        &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
    };
    const int nf = (int)(sizeof(fields) / sizeof(fields[0]));

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                double lo = c43*F.x[i][j][1]    - c13*F.x[i][j][2];
                double hi = c43*F.x[i][j][km-2] - c13*F.x[i][j][km-3];
                F.x[i][j][0] = F.x[i][j][km-1] = 0.5*(lo + hi);
            }
        }
    }

    // --- (4) Shapiro damping across the phi seam. ---
    // k=0 and k=km-1 are the SAME physical longitude and are not evolved by the RK4 (its phi
    // loop runs k=1..km-2); the reconstruction above pins them to 0.5*(x[1]+x[km-2]). Because
    // that slaves the seam value to its own neighbours, the discrete d2/dphi2 self-damping at
    // the seam-adjacent cells k=1 and k=km-2 drops from -2 to -1.5 — a 25% loss of numerical
    // zonal diffusion exactly at the seam. Combined with the 1/sin^2(theta) metric that leaves
    // an under-damped zonal mode that runs away in ATOM. In ATJUP no such mode is present yet
    // (see the seam_damp() note above), which is why this defaults to 0 passes.
    //
    // One explicit 1-2-1 pass with coefficient 0.25 fully removes the 2*dphi mode in steady
    // state; two passes extend the reach to ~4*dphi, which is needed when a jet crossing the
    // seam regenerates seam energy faster than a single pass absorbs it. Only u,v,w are damped
    // (restoreVar copies them to un,vn,wn after all BCs, so the n-level follows automatically),
    // and solid SeaMount cells act as no-flux: a solid neighbour contributes the cell's own
    // value, relaxing it toward the seam rather than dragging it toward zero.
    const int npass = BCJupKnobs::seam_damp();
    if(npass > 0){
        constexpr double seam_coeff = 0.25;
        const int km2 = km - 2;
        const int km3 = km - 3;

        #pragma omp parallel for schedule(static)
        for(int i = 0; i < im; i++){
            for(int j = 0; j < jm; j++){
                auto is_fluid = [&](int kk){ return m.SeaMount.x[i][j][kk] != 1.0; };

                auto smooth_seam = [&](Array& F){
                    const double f_km2 = F.x[i][j][km2];
                    const double f_0   = F.x[i][j][0];
                    const double f_1   = F.x[i][j][1];

                    const bool fl_km2 = is_fluid(km2);
                    const bool fl_0   = is_fluid(0);
                    const bool fl_1   = is_fluid(1);

                    // Periodic neighbours; substitute the cell's own value at a solid neighbour.
                    const double w_km2 = is_fluid(km3) ? F.x[i][j][km3] : f_km2;
                    const double e_km2 = fl_0 ? f_0   : f_km2;   // east of km-2 is the seam
                    const double w_0   = fl_km2 ? f_km2 : f_0;   // west of the seam is km-2
                    const double e_0   = fl_1 ? f_1   : f_0;     // east of the seam is k=1
                    const double w_1   = fl_0 ? f_0   : f_1;     // west of k=1 is the seam
                    const double e_1   = is_fluid(2) ? F.x[i][j][2] : f_1;

                    if(fl_km2)
                        F.x[i][j][km2] = f_km2 + seam_coeff*(w_km2 - 2.0*f_km2 + e_km2);
                    if(fl_0)
                        F.x[i][j][0] = F.x[i][j][km-1] =
                            f_0 + seam_coeff*(w_0 - 2.0*f_0 + e_0);
                    if(fl_1)
                        F.x[i][j][1] = f_1 + seam_coeff*(w_1 - 2.0*f_1 + e_1);
                };

                for(int p = 0; p < npass; p++){
                    smooth_seam(m.u);
                    smooth_seam(m.v);
                    smooth_seam(m.w);
                }
            }
        }
    }
}


inline void BC_Jup::bcSeaMount()
{
    cout << endl << "      ATJUP: BC_seamount" << endl;

    const int im = m.im;
    const int jm = m.jm;
    const int km = m.km;

    int j_ellipse = 0;
    int a = im-1;
    int b = im-1;
    const int j_0 = 112;  // center of Great Red Spot, 22° south of equator
    const int k_0 = 180;
    const int i_0 = 35;

    // DIAGNOSTIC KNOB, default off: build no obstacle at all, leaving a smooth spherical
    // shell. The secular growth of max|w| is anchored to a single fluid cell two cells
    // outside the staircase flank of this cone, so a run without the cone is the control
    // that separates "the obstacle edge makes it" from "the momentum budget makes it".
    // i_topography is still filled below, all zeros, so every consumer stays valid.
    const bool no_mount = BCJupKnobs::env_int("ATJUP_NO_SEAMOUNT", 0) != 0;

    for(int i = 0; !no_mount && i < i_0; i++){
        if(i <= 10) a = b = im-1;
        for(int k = k_0-a; k <= k_0+a; k++){
            for(int j = j_0-b; j <= j_0+b; j++){
                if(j <= j_0){
                    j_ellipse = j_0 - (int)sqrt(pow(b,2) - pow((b*(k-k_0)/a),2));
                    if(j >= j_ellipse)  m.SeaMount.x[i][j][k] = 1.0;
                } else {
                    j_ellipse = j_0 + (int)sqrt(pow(b,2) - pow((b*(k-k_0)/a),2));
                    if(j <= j_ellipse)  m.SeaMount.x[i][j][k] = 1.0;
                }
                if(i == i_0) m.Topography.y[j][k] = m.SeaMount.x[i_0][j][k];
            }
        }
        a = (im-1)-i;
        b = (im-1)-i;
    }

    m.i_topography.assign(jm, std::vector<int>(km, 0));
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_surf = 0;
            while(i_surf < im && m.SeaMount.x[i_surf][j][k] == 1.0) ++i_surf;
            m.i_topography[j][k] = i_surf;
        }
    }

    cout << "      ATJUP: BC_seamount ended" << endl;
}


inline void BC_Jup::bcSolidGround()
{
    cout << endl << "      ATJUP: BC_solidground" << endl;

    const int im = m.im, jm = m.jm, km = m.km;

    // All fields receive the same treatment: zero at interior solid cells,
    // cubic extrapolation from the fluid side at surface solid cells.
    // bcSolidGround is called after RungeKuttaJup() completes — diagnostic
    // fields (massflux, difflux, forces, Q) are frozen throughout the RK4
    // substeps and are in a fully computed state here, so extrapolating them
    // is as valid as extrapolating the primary state variables.
    Array* scalars[] = {
        &m.t, &m.p_stat, &m.p_dyn,
        &m.ch4, &m.ch4_cloud, &m.ch4_ice,
        &m.h2o, &m.h2o_cloud, &m.h2o_ice,
        &m.h2s,
        &m.nh3, &m.nh3_cloud, &m.nh3_ice,
        &m.nh4sh,
        &m.j_h2s,  &m.j_nh3,  &m.j_nh4sh,
        &m.jT_h2s, &m.jT_nh3, &m.jT_nh4sh,
        &m.w_h2s,  &m.w_nh3,  &m.w_nh4sh,
        &m.massflux_h2s,  &m.massflux_nh3,  &m.massflux_nh4sh,
        &m.fluxlim_nh4sh,
        &m.difflux_h2s,   &m.difflux_nh3,   &m.difflux_nh4sh,
        &m.thermalmassflux,
        &m.CoriolisForce, &m.CentrifugalForce,
        &m.BuoyancyForce, &m.PresGradForce,
        &m.Q_Latent, &m.Q_Sensible
    };
    const int ns = (int)(sizeof(scalars) / sizeof(scalars[0]));
    const double c43 = m.c43, c13 = m.c13;

    // t and p_stat (entries 0 and 1 above) are THERMODYNAMIC STATE, not a transported
    // quantity that can be "absent" in a solid. Zeroing them in the obstacle interior, as
    // this loop used to do for every field, put T = 0 K and p = 0 bar into those cells, and
    // every routine that forms a thermodynamic expression divides by one of them:
    //   Forces()/rhs_u buoyancy   g*(p_stat+p_dyn) / (r_mix*R_mix*t*t_ref)  -> x/0 = inf
    //   ChemistryJup::DiffMassFluxJup  LT_x / t                             -> x/0 = inf
    //   SaturationAdjustmentJup   exp(-L0/T + ...) and E/p_u                -> inf, NaN
    //   Thermo_Jup::Latent_Heat   E/(p_u - E)                               -> NaN
    // A 0 K cell also puts a ~120 K artificial jump into the temperature Laplacian of any
    // fluid cell touching the obstacle face, which is a physics error independent of the
    // NaN. The interior of the obstacle never enters the flow solution (RungeKuttaJup skips
    // SeaMount cells), so leaving t and p_stat at the ambient hydrostatic profile they were
    // initialised with is both harmless and physically the right filler. Everything from
    // p_dyn onwards is genuinely zero inside a solid and is still zeroed below.
    const int n_state = 2;     // scalars[0] = t, scalars[1] = p_stat

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(m.SeaMount.x[i][j][k] != 1.0) continue;

                // u, v, w are zero everywhere inside the SeaMount and on its surface.
                m.u.x[i][j][k] = 0.0;
                m.v.x[i][j][k] = 0.0;
                m.w.x[i][j][k] = 0.0;

                // Lambdas defined inside the loop body so GCC can track their use
                // through the OpenMP outlined function (avoids -Wunused-but-set-variable).

                // Out-of-bounds counts as solid (domain walls are not open air).
                auto solid = [&](int ii, int jj, int kk) -> bool {
                    if(ii < 0 || ii >= im || jj < 0 || jj >= jm || kk < 0 || kk >= km) return true;
                    return m.SeaMount.x[ii][jj][kk] == 1.0;
                };
                // True only when (ii,jj,kk) is in-bounds and fluid.
                auto fluid = [&](int ii, int jj, int kk) -> bool {
                    if(ii < 0 || ii >= im || jj < 0 || jj >= jm || kk < 0 || kk >= km) return false;
                    return m.SeaMount.x[ii][jj][kk] != 1.0;
                };
                // 2-point extrapolation onto (i,j,k): f[s] = (4/3)f[a] - (1/3)f[b].
                // Consistent with bcScalarSurfSur and bcPhi. The 3-point cubic
                // (3f[a]-3f[b]+f[c]) amplifies alternating errors by up to 7x per step
                // and causes blow-up near the obstacle surface around iter_n=32.
                // The 2-point formula limits amplification to (5/3)x, which the fluid
                // diffusion keeps bounded.
                // ATJUP_BC_MOUNT_COPY=1 replaces it by the plain copy f[s] = f[a], which cannot
                // overshoot — see the knob's comment for the flank-cooling measurement it is
                // meant to settle.
                const bool mount_copy = BCJupKnobs::mount_copy() != 0;
                auto extrap = [&](int ia, int ja, int ka,
                                  int ib, int jb, int kb) {
                    for(int f = 0; f < ns; f++){
                        double*** x = scalars[f]->x;
                        x[i][j][k] = mount_copy
                                   ? x[ia][ja][ka]
                                   : c43*x[ia][ja][ka] - c13*x[ib][jb][kb];
                    }
                };

                // Surface cell: at least one face-neighbor is fluid.
                // Interior cell: all six face-neighbors are solid or at a domain wall.
                const bool is_surface =
                    !solid(i+1,j,k) || !solid(i-1,j,k) ||
                    !solid(i,j+1,k) || !solid(i,j-1,k) ||
                    !solid(i,j,k+1) || !solid(i,j,k-1);

                if(!is_surface){
                    // from n_state on: t and p_stat keep their ambient values (see above)
                    for(int f = n_state; f < ns; f++)
                        scalars[f]->x[i][j][k] = 0.0;
                    continue;
                }

                // Surface: 2-point extrapolation in the first direction with two
                // consecutive in-bounds fluid cells. Priority: ±i, ±j, ±k.
                // Fallback to zero when no such pair exists (thin spur cells).
                if(fluid(i+1,j,k) && fluid(i+2,j,k)){
                    extrap(i+1,j,k, i+2,j,k);
                } else if(fluid(i-1,j,k) && fluid(i-2,j,k)){
                    extrap(i-1,j,k, i-2,j,k);
                } else if(fluid(i,j+1,k) && fluid(i,j+2,k)){
                    extrap(i,j+1,k, i,j+2,k);
                } else if(fluid(i,j-1,k) && fluid(i,j-2,k)){
                    extrap(i,j-1,k, i,j-2,k);
                } else if(fluid(i,j,k+1) && fluid(i,j,k+2)){
                    extrap(i,j,k+1, i,j,k+2);
                } else if(fluid(i,j,k-1) && fluid(i,j,k-2)){
                    extrap(i,j,k-1, i,j,k-2);
                } else {
                    // Thin spur cell with no usable fluid pair to extrapolate from: same
                    // treatment as the interior, t and p_stat keep their ambient values.
                    for(int f = n_state; f < ns; f++)
                        scalars[f]->x[i][j][k] = 0.0;
                }
            }
        }
    }

    // ---- Final pass: bit-level sanitisation of the turbulence fields ----
    // Ported from ATOM's BC_Atm.h "Pass 5". Now that k* and dis* are prognostic
    // (rhs_tke / rhs_dis), a single non-finite cell anywhere in the closure propagates:
    // TurbulenceJup::compute_vel_star reads v,w at the first fluid layer, and a non-finite
    // velocity there makes u_tau infinite, which the ABL background floor k_bg then writes
    // straight into k*. Observed at 87 deg N, i=1, iteration 2. Clip every turbulence array
    // to a physically generous ceiling and reset non-finite values to a sane default, as the
    // last operation before the next RK4 cycle reads them.
    //
    // Bit-level test rather than std::isfinite/std::max: with -ffast-math the compiler may
    // assume the operands are finite, and ATOM measured k reaching 7e200 despite an RK4 cap
    // for exactly that reason. Reading the IEEE-754 exponent field has no such hazard.
    //
    // This touches ONLY the turbulence arrays, which are identically zero unless ATJUP_TURB
    // is set, so a run without the closure is bit-identical. The non-finite VELOCITY cell that
    // triggers this at high latitude is a separate, pre-existing ATJUP problem and is
    // deliberately not masked here.
    {
        const double tke_max_nd = 1000.0 / (m.u_0 * m.u_0);        // 1000 m2/s2, as in RK4
        const double nue_max    = 1.0e5 / (m.u_0 * m.L_atm * 1.0e3); // ATJUP_NUE_MAX default
        const double prod_max   = 1.0e4;                            // generous; real prod ~ 1
        const double dis_max    = 1.0e6;                            // generous for omega*

        auto safe_clamp = [](double v, double lo, double hi) -> double {
            std::uint64_t bits;
            std::memcpy(&bits, &v, sizeof(bits));
            if((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) return lo;
            return (v < lo) ? lo : ((v > hi) ? hi : v);
        };

        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < im; i++){
            for(int j = 0; j < jm; j++){
                for(int k = 0; k < km; k++){
                    m.tke.x[i][j][k]        = safe_clamp(m.tke.x[i][j][k],        0.0,       tke_max_nd);
                    m.tken.x[i][j][k]       = safe_clamp(m.tken.x[i][j][k],       0.0,       tke_max_nd);
                    m.dis.x[i][j][k]        = safe_clamp(m.dis.x[i][j][k],        1.0e-10,   dis_max);
                    m.disn.x[i][j][k]       = safe_clamp(m.disn.x[i][j][k],       1.0e-10,   dis_max);
                    m.nue.x[i][j][k]        = safe_clamp(m.nue.x[i][j][k],        0.0,       nue_max);
                    m.prod.x[i][j][k]       = safe_clamp(m.prod.x[i][j][k],       0.0,       prod_max);
                    m.tke_source.x[i][j][k] = safe_clamp(m.tke_source.x[i][j][k], -prod_max, prod_max);
                    m.dis_source.x[i][j][k] = safe_clamp(m.dis_source.x[i][j][k], -prod_max, prod_max);
                }
            }
        }
    }

    cout << "      ATJUP: BC_solidground ended" << endl;
}


inline void BC_Jup::bcVelSurfSur()
{
    cout << endl << "      ATJUP: BC_vel_surf_sur" << endl;

    const int im = m.im, jm = m.jm, km = m.km;

    constexpr double coeff  = 0.01;
    constexpr double coeff5 = 0.9;

    auto is_land = [&](int i, int j, int k) { return m.SeaMount.x[i][j][k] == 1.0; };
    auto is_air  = [&](int i, int j, int k) { return m.SeaMount.x[i][j][k] != 1.0; };

    #pragma omp parallel for schedule(static)
    for (int i = 1; i < im-1; i++) {
        for (int j = 1; j < jm-1; j++) {
            for (int k = 1; k < km-1; k++) {

                if (!is_land(i, j, k)) continue;

                // i-direction (radial: SeaMount → air outward)
                if (i < im-2 && is_air(i+1, j, k)) {
                    m.u.x[i][j][k]   = 0.0;
                    m.u.x[i+1][j][k] = 0.0;

                    m.v.x[i][j][k]    = 0.0;
                    m.v.x[i+1][j][k] *= coeff;
                    m.v.x[i+2][j][k] *= coeff5;

                    m.w.x[i][j][k]    = 0.0;
                    m.w.x[i+1][j][k] *= coeff;
                    m.w.x[i+2][j][k] *= coeff5;
                }

                // j-direction (meridional: SeaMount → air south/north)
                if (j < jm-2 && is_air(i, j+1, k) && is_air(i, j+2, k)) {
                    m.u.x[i][j][k]    = 0.0;
                    m.u.x[i][j+1][k] *= coeff;
                    m.u.x[i][j+2][k] *= coeff5;

                    m.v.x[i][j][k]    = 0.0;
                    m.v.x[i][j+1][k]  = 0.0;

                    m.w.x[i][j][k]    = 0.0;
                    m.w.x[i][j+1][k] *= coeff;
                    m.w.x[i][j+2][k] *= coeff5;
                }
                if (j >= 2 && is_air(i, j-1, k) && is_air(i, j-2, k)) {
                    m.u.x[i][j][k]    = 0.0;
                    m.u.x[i][j-1][k] *= coeff;
                    m.u.x[i][j-2][k] *= coeff5;

                    m.v.x[i][j][k]    = 0.0;
                    m.v.x[i][j-1][k]  = 0.0;

                    m.w.x[i][j][k]    = 0.0;
                    m.w.x[i][j-1][k] *= coeff;
                    m.w.x[i][j-2][k] *= coeff5;
                }

                // k-direction (zonal: SeaMount → air east/west)
                if (k < km-2 && is_air(i, j, k+1) && is_air(i, j, k+2)) {
                    m.u.x[i][j][k]    = 0.0;
                    m.u.x[i][j][k+1] *= coeff;
                    m.u.x[i][j][k+2] *= coeff5;

                    m.v.x[i][j][k]    = 0.0;
                    m.v.x[i][j][k+1] *= coeff;
                    m.v.x[i][j][k+2] *= coeff5;

                    m.w.x[i][j][k]    = 0.0;
                    m.w.x[i][j][k+1]  = 0.0;
                }
                if (k >= 2 && is_air(i, j, k-1) && is_air(i, j, k-2)) {
                    m.u.x[i][j][k]    = 0.0;
                    m.u.x[i][j][k-1] *= coeff;
                    m.u.x[i][j][k-2] *= coeff5;

                    m.v.x[i][j][k]    = 0.0;
                    m.v.x[i][j][k-1] *= coeff;
                    m.v.x[i][j][k-2] *= coeff5;

                    m.w.x[i][j][k]    = 0.0;
                    m.w.x[i][j][k-1]  = 0.0;
                }
            }
        }
    }
    cout << "      ATJUP: BC_vel_surf_sur ended" << endl;
}


inline void BC_Jup::bcScalarSurfSur()
{
    cout << endl << "      ATJUP: BC_scalar_surf_sur" << endl;

    const int im = m.im, jm = m.jm, km = m.km;
    const double c43 = m.c43, c13 = m.c13;

    auto is_land = [&](int i, int j, int k) { return m.SeaMount.x[i][j][k] == 1.0; };
    auto is_air  = [&](int i, int j, int k) { return m.SeaMount.x[i][j][k] != 1.0; };

    Array* scalars[] = {
        &m.u, &m.v, &m.w,
        &m.t, &m.p_stat,
        &m.ch4, &m.ch4_cloud, &m.ch4_ice,
        &m.h2o, &m.h2o_cloud, &m.h2o_ice,
        &m.h2s, &m.j_h2s, &m.jT_h2s,
        &m.nh3, &m.nh3_cloud, &m.nh3_ice,
        &m.j_nh3, &m.jT_nh3,
        &m.nh4sh, &m.j_nh4sh, &m.jT_nh4sh,
        &m.difflux_h2s, &m.difflux_nh3, &m.difflux_nh4sh, 
        &m.massflux_h2s, &m.massflux_nh3, &m.massflux_nh4sh, 
        &m.w_h2s, &m.w_nh3, &m.w_nh4sh, 
        &m.thermalmassflux,
        &m.Q_Latent, &m.Q_Sensible,
        &m.BuoyancyForce, &m.CoriolisForce, &m.CentrifugalForce, &m.PresGradForce,
    };
    constexpr int ns = (int)(sizeof(scalars) / sizeof(scalars[0]));

    // Neumann (zero normal gradient) BC at SeaMount–air interfaces.
    // Vertical (i) direction is preferred at corners.
    #pragma omp parallel for schedule(static)
    for (int i = 1; i < im-1; i++) {
        for (int j = 1; j < jm-1; j++) {
            for (int k = 1; k < km-1; k++) {

                if (!is_land(i, j, k)) continue;

                // ---- i-direction (preferred at corners) ----
                if (is_air(i+1, j, k)) {
                    if (i + 2 < im) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = c43 * x[i+1][j][k] - c13 * x[i+2][j][k];
                        }
                    } else {
                        for (int f = 0; f < ns; f++)
                            scalars[f]->x[i][j][k] = scalars[f]->x[i+1][j][k];
                    }
                    continue;
                }

                // ---- j-direction ----
                if (j < jm-2 && is_air(i, j+1, k) && is_air(i, j+2, k)) {
                    for (int f = 0; f < ns; f++) {
                        double*** x = scalars[f]->x;
                        x[i][j][k] = c43 * x[i][j+1][k] - c13 * x[i][j+2][k];
                    }
                } else if (j >= 2 && is_air(i, j-1, k) && is_air(i, j-2, k)) {
                    for (int f = 0; f < ns; f++) {
                        double*** x = scalars[f]->x;
                        x[i][j][k] = c43 * x[i][j-1][k] - c13 * x[i][j-2][k];
                    }
                }

                // ---- k-direction ----
                if (k < km-2 && is_air(i, j, k+1) && is_air(i, j, k+2)) {
                    for (int f = 0; f < ns; f++) {
                        double*** x = scalars[f]->x;
                        x[i][j][k] = c43 * x[i][j][k+1] - c13 * x[i][j][k+2];
                    }
                } else if (k >= 2 && is_air(i, j, k-1) && is_air(i, j, k-2)) {
                    for (int f = 0; f < ns; f++) {
                        double*** x = scalars[f]->x;
                        x[i][j][k] = c43 * x[i][j][k-1] - c13 * x[i][j][k-2];
                    }
                }
            }
        }
    }

    // Project SeaMount surface values down to i=0 so the base layer always
    // carries the topmost fluid-cell values (equivalent to ATOM's i_topography).
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            // Find first air cell from the bottom up.
            int i_surf = 0;
            while (i_surf < im && is_land(i_surf, j, k)) ++i_surf;
            if (i_surf == 0) continue;  // no SeaMount at this (j,k)
            if (i_surf >= im) continue; // fully blocked column (degenerate)

            m.t.x[0][j][k]         = m.t.x[i_surf][j][k];
            m.p_stat.x[0][j][k]    = m.p_stat.x[i_surf][j][k];
            m.ch4.x[0][j][k]       = m.ch4.x[i_surf][j][k];
            m.ch4_cloud.x[0][j][k] = m.ch4_cloud.x[i_surf][j][k];
            m.ch4_ice.x[0][j][k]   = m.ch4_ice.x[i_surf][j][k];
            m.h2o.x[0][j][k]       = m.h2o.x[i_surf][j][k];
            m.h2o_cloud.x[0][j][k] = m.h2o_cloud.x[i_surf][j][k];
            m.h2o_ice.x[0][j][k]   = m.h2o_ice.x[i_surf][j][k];
            m.h2s.x[0][j][k]       = m.h2s.x[i_surf][j][k];
            m.nh3.x[0][j][k]       = m.nh3.x[i_surf][j][k];
            m.nh3_cloud.x[0][j][k] = m.nh3_cloud.x[i_surf][j][k];
            m.nh3_ice.x[0][j][k]   = m.nh3_ice.x[i_surf][j][k];
            m.nh4sh.x[0][j][k]     = m.nh4sh.x[i_surf][j][k];
        }
    }

    cout << "      ATJUP: BC_scalar_surf_sur ended" << endl;
}


inline void BC_Jup::initTropopauseLayers()
{
    m.tropopause_layers = std::vector<double>(m.jm, m.tropopause_pole);
    cout << endl << "      ATJUP: init_tropopause_layers" << endl;

    const int i_max  = m.im - 1;
    const int j_max  = m.jm - 1;
    const int j_half = j_max / 2;
    const double coeff_pole = 285.0;

    for(int j = j_half; j >= 0; j--){
        double x = coeff_pole * (1.0 - (double)(j_half - j) / (double)j_half);
        m.tropopause_layers[j] = JupiterUtils::Agnesi(m.tropopause_equator, x);
        m.tropopause_layers[j] = std::round(m.tropopause_layers[j]
            / m.L_atm * (double)i_max);
        m.tropopause_layers[j] = m.tropopause_equator / m.L_atm * (double)i_max;
    }

    for(int j = j_max; j > j_half; j--)
        m.tropopause_layers[j] = m.tropopause_layers[j_max - j];

    cout << "      ATJUP: init_tropopause_layers ended" << endl;
}
