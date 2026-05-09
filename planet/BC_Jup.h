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
        &m.Q_Latent, &m.Q_Sensible
    };
    const int nf = (int)(sizeof(fields) / sizeof(fields[0]));

    // 2-point Neumann extrapolation: f[s] = (4/3)f[a] - (1/3)f[b].
    // The 3-point cubic (3f[a]-3f[b]+f[c]) amplifies alternating errors by 7x
    // per call and blows up near the SeaMount contour (same reason bcSolidGround
    // was switched to the 2-point formula; bcRadius had the same latent bug).
    #pragma omp parallel for schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                F.x[0][j][k]    = c43*F.x[1][j][k]    - c13*F.x[2][j][k];
                F.x[im-1][j][k] = c43*F.x[im-2][j][k] - c13*F.x[im-3][j][k];
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
        &m.Q_Latent, &m.Q_Sensible
    };
    const int nf = (int)(sizeof(extrap_fields) / sizeof(extrap_fields[0]));

    // Flux fields that are singular-prone near poles: zero at pole boundary.
    Array* zero_at_poles[] = {
        &m.massflux_h2s, &m.massflux_nh3, &m.massflux_nh4sh,
        &m.fluxlim_nh4sh,
        &m.difflux_h2s,  &m.difflux_nh3,  &m.difflux_nh4sh,
    };
    const int nz = (int)(sizeof(zero_at_poles) / sizeof(zero_at_poles[0]));

    #pragma omp parallel for schedule(static)
//    for(int k = 1; k < km-1; k++){
//        for(int i = 1; i < im-1; i++){
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            m.v.x[i][0][k]    = 0.0;
            m.v.x[i][jm-1][k] = 0.0;
            m.w.x[i][0][k]    = 0.0;
            m.w.x[i][jm-1][k] = 0.0;

            for(int f = 0; f < nf; f++){
                Array& F = *extrap_fields[f];
                F.x[i][0][k]    = c43*F.x[i][1][k]    - c13*F.x[i][2][k];
                F.x[i][jm-1][k] = c43*F.x[i][jm-2][k] - c13*F.x[i][jm-3][k];
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
        &m.Q_Latent, &m.Q_Sensible
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

    for(int i = 0; i < i_0; i++){
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
                auto extrap = [&](int ia, int ja, int ka,
                                  int ib, int jb, int kb) {
                    for(int f = 0; f < ns; f++){
                        double*** x = scalars[f]->x;
                        x[i][j][k] = c43*x[ia][ja][ka] - c13*x[ib][jb][kb];
                    }
                };

                // Surface cell: at least one face-neighbor is fluid.
                // Interior cell: all six face-neighbors are solid or at a domain wall.
                const bool is_surface =
                    !solid(i+1,j,k) || !solid(i-1,j,k) ||
                    !solid(i,j+1,k) || !solid(i,j-1,k) ||
                    !solid(i,j,k+1) || !solid(i,j,k-1);

                if(!is_surface){
                    for(int f = 0; f < ns; f++)
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
                    for(int f = 0; f < ns; f++)
                        scalars[f]->x[i][j][k] = 0.0;
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
