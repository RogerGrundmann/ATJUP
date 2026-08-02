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
#include "BoundaryConditions.h"
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

    // The five DOMAIN boundary knobs that used to live here — rigid_lid, top_taper,
    // pole_copy, radius_copy, seam_damp and the lid temperature pin — moved into the
    // shared BoundaryConditions.h along with the routines they gate, and their defaults
    // are now cJupiterModel::bc_default_*(). Their names and values are unchanged:
    // ATJUP_BC_RIGID_LID and ATJUP_BC_TOP_TAPER still default to 1, ATJUP_BC_POLE_COPY
    // to 1, ATJUP_BC_RADIUS_COPY and ATJUP_BC_SEAM_DAMP to 0.
    //
    // What stays is the OBSTACLE knob, because the obstacle routines below are ATJUP's
    // alone — ATSAT has no solid body at all.

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
}


// ---------------------------------------------------------------------------
// The three DOMAIN boundary conditions are the shared BoundaryConditions<Planet> template.
// What ATJUP supplies is WHICH FIELDS each one treats; everything else it needs is answered by
// cJupiterModel::bc_margin(), bc_default_form() and bc_default_*().
// ---------------------------------------------------------------------------
inline void BC_Jup::bcRadius(){ BoundaryConditions<cJupiterModel>(m).bcRadius(); }
inline void BC_Jup::bcTheta() { BoundaryConditions<cJupiterModel>(m).bcTheta();  }
inline void BC_Jup::bcPhi()   { BoundaryConditions<cJupiterModel>(m).bcPhi();    }


// k* and dis* are prognostic (RHS_Jup_Turb.cpp), so their domain boundaries must be extrapolated
// like every other transported field. Without this they stay frozen at the array initial values
// (k*=0, dis*=1e-10) while the interior runs at dis* ~ 16, and the resulting permanent Laplacian
// at i=0/im-1, the poles and the phi seam drains dis* to its floor within ~20 iterations; nue*
// then saturates its ceiling and, with ATJUP_TURB_COUPLING on, the eddy viscosity destroys the
// momentum field. ATJUP carries them in the main lists below. ATSAT instead gives them their own
// clamped, 2-point treatment, which is the better of the two — see BoundaryConditions.h.
inline std::vector<Array*> cJupiterModel::bc_fields_radius(){
    return {
        &t, &u, &v, &w,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,  &j_nh3,  &j_nh4sh,
        &jT_h2s, &jT_nh3, &jT_nh4sh,
        &w_h2s,  &w_nh3,  &w_nh4sh,
        &massflux_h2s,  &massflux_nh3,  &massflux_nh4sh,
        &difflux_h2s,   &difflux_nh3,   &difflux_nh4sh,
        &thermalmassflux,
        &CoriolisForce, &CentrifugalForce,
        &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible,
        &tke, &dis, &nue, &prod, &tke_source, &dis_source
    };
}

// All fields except v, w and the flux fields receive the pole extrapolation.
inline std::vector<Array*> cJupiterModel::bc_fields_theta_extrap(){
    return {
        &t, &u,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,  &j_nh3,  &j_nh4sh,
        &jT_h2s, &jT_nh3, &jT_nh4sh,
        &w_h2s,  &w_nh3,  &w_nh4sh,
        &thermalmassflux,
        &CoriolisForce, &CentrifugalForce,
        &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible,
        &tke, &dis, &nue, &prod, &tke_source, &dis_source
    };
}

// v and w have no meaning on the axis; the mass and diffusive fluxes are singular-prone there.
// NOTE ATSAT EXTRAPOLATES massflux_* and difflux_* at the poles instead of zeroing them. That is
// a real disagreement about the physics, it is preserved rather than resolved by sharing, and it
// is visible now only because both models finally state their lists in one place each.
inline std::vector<Array*> cJupiterModel::bc_fields_theta_zero(){
    return {
        &v, &w,
        &massflux_h2s, &massflux_nh3, &massflux_nh4sh,
        &fluxlim_nh4sh,
        &difflux_h2s,  &difflux_nh3,  &difflux_nh4sh
    };
}

// NOTE this list omits j_nh4sh and jT_nh4sh, which ATSAT's phi list carries. Preserved as found.
inline std::vector<Array*> cJupiterModel::bc_fields_phi(){
    return {
        &t, &u, &v, &w,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,  &j_nh3,
        &jT_h2s, &jT_nh3,
        &w_h2s,  &w_nh3,  &w_nh4sh,
        &massflux_h2s,  &massflux_nh3,  &massflux_nh4sh,
        &fluxlim_nh4sh,
        &difflux_h2s,   &difflux_nh3,   &difflux_nh4sh,
        &thermalmassflux,
        &CoriolisForce, &CentrifugalForce,
        &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible,
        &tke, &dis, &nue, &prod, &tke_source, &dis_source
    };
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
