#pragma once

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Array.h"
#include "Array_1D.h"

// Full definition of cJupiterModel is provided by the deferred include at the
// bottom of cJupiterModel.h.  That file forward-declares ChemistryJup before
// defining cJupiterModel, then includes this header afterwards so the inline
// bodies below can reference m.xxx freely.
#include "cJupiterModel.h"

class ChemistryJup {
    friend class cJupiterModel;

public:
    explicit ChemistryJup(cJupiterModel& model) : m(model) {}
    ~ChemistryJup() = default;
    ChemistryJup(const ChemistryJup&) = delete;
    ChemistryJup& operator=(const ChemistryJup&) = delete;



    void ChemMassRateJup()
    {
        using namespace std;
        cout << endl << "      ATJUP: ChemMassRateJup" << endl;

        auto begin = chrono::high_resolution_clock::now();

        const double keq = m.m_nh4sh / (m.m_nh3 * m.m_h2s);
        const double A   = 15000.0;
        const double B   = -0.5;
        const double T_d = 3020.0;

        // Zero the polar ghost bands (j < 3 and j > jm-4) that lie outside the
        // RungeKutta domain (which runs j = 3..jm-4).  Concentrations there are
        // only touched by SaturationAdjust, so the quadratic reaction term
        // kf·c_NH3·c_H2S has no dynamical sink and grows without bound.
        // The massflux values produced there are never consumed by the RHS, so
        // they show up only as a spurious growing artifact in output.
        #ifdef _OPENMP
        #pragma omp parallel for collapse(2) schedule(static)
        #endif
        for (int i = 0; i < m.im; i++) {
            for (int k = 0; k < m.km; k++) {
                for (int j = 0; j < 3; j++) {
                    m.w_nh3.x[i][j][k]   = m.w_h2s.x[i][j][k]   = m.w_nh4sh.x[i][j][k]   = 0.0;
                    m.massflux_nh3.x[i][j][k] = m.massflux_h2s.x[i][j][k] = m.massflux_nh4sh.x[i][j][k] = 0.0;
                }
                for (int j = m.jm-3; j < m.jm; j++) {
                    m.w_nh3.x[i][j][k]   = m.w_h2s.x[i][j][k]   = m.w_nh4sh.x[i][j][k]   = 0.0;
                    m.massflux_nh3.x[i][j][k] = m.massflux_h2s.x[i][j][k] = m.massflux_nh4sh.x[i][j][k] = 0.0;
                }
            }
        }

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (int k = 1; k < m.km-1; k++) {
            for (int j = 3; j < m.jm-3; j++) {
                for (int i = 1; i < m.im-1; i++) {

                    // Interior solid cells (no fluid face-neighbor): zero and skip.
                    // Surface solid cells (at least one fluid face-neighbor) keep their
                    // extrapolated values so adjacent fluid cells see correct gradients.
                    if (m.SeaMount.x[i][j][k] == 1.0
                        && m.SeaMount.x[i-1][j][k] == 1.0 && m.SeaMount.x[i+1][j][k] == 1.0
                        && m.SeaMount.x[i][j-1][k] == 1.0 && m.SeaMount.x[i][j+1][k] == 1.0
                        && m.SeaMount.x[i][j][k-1] == 1.0 && m.SeaMount.x[i][j][k+1] == 1.0) {
                        m.w_nh3.x[i][j][k]          = 0.0;
                        m.w_h2s.x[i][j][k]          = 0.0;
                        m.w_nh4sh.x[i][j][k]        = 0.0;
                        m.massflux_nh3.x[i][j][k]   = 0.0;
                        m.massflux_h2s.x[i][j][k]   = 0.0;
                        m.massflux_nh4sh.x[i][j][k] = 0.0;
                        continue;
                    }

                    const double t_u = m.t.x[i][j][k] * m.t_ref;
                    const double kf  = react_rate_const(t_u, T_d, A, B);
                    const double kb  = kf / keq;

                    if ((t_u <= m.t_0_nh4sh) && (t_u >= m.t_00_nh4sh)) {
                        const double sum_c   = m.nh3.x[i][j][k]   / m.m_nh3
                                             + m.h2s.x[i][j][k]   / m.m_h2s
                                             + m.nh4sh.x[i][j][k] / m.m_nh4sh;
                        const double denom   = (sum_c == 0.0) ? 0.0 : m.r_mix / sum_c;

                        const double c_nh3   = m.nh3.x[i][j][k]   / m.m_nh3   * denom;
                        const double c_h2s   = m.h2s.x[i][j][k]   / m.m_h2s   * denom;
                        const double c_nh4sh = m.nh4sh.x[i][j][k] / m.m_nh4sh * denom;

                        const double R_diff = kf * c_nh3 * c_h2s - kb * c_nh4sh;

                        m.w_nh3.x[i][j][k]   = -m.m_nh3   * R_diff;
                        m.w_h2s.x[i][j][k]   = -m.m_h2s   * R_diff;
                        m.w_nh4sh.x[i][j][k] =  m.m_nh4sh * R_diff;
                    }

                    m.massflux_h2s.x[i][j][k]   = m.w_h2s.x[i][j][k]   - m.difflux_h2s.x[i][j][k];
                    m.massflux_nh3.x[i][j][k]   = m.w_nh3.x[i][j][k]   - m.difflux_nh3.x[i][j][k];
                    m.massflux_nh4sh.x[i][j][k] = m.w_nh4sh.x[i][j][k] - m.difflux_nh4sh.x[i][j][k];
                }
            }
        }

        auto end     = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ChemMassRateJup\n", elapsed.count() * 1e-9);
        cout << "      ATJUP: ChemMassRateJup ended" << endl;
    }

    void DiffMassFluxJup()
    {
        using namespace std;
        cout << endl << "      ATJUP: DiffMassFluxJup" << endl;

        auto begin = chrono::high_resolution_clock::now();

        const double DT_nh3   = m.mue_nh3   / (m.rg_nh3   * m.sc_nh3);
        const double DT_h2s   = m.mue_h2s   / (m.rg_h2s   * m.sc_h2s);
        const double DT_nh4sh = m.mue_nh4sh / (m.rg_nh4sh * m.sc_nh4sh);
        const double D_nh3    = DT_nh3;
        const double D_h2s    = DT_h2s;
        const double D_nh4sh  = DT_nh4sh;

        const double M_mix    = m.X_h2 * m.m_h2 + m.X_he * m.m_he
                              + m.X_h2o * m.m_h2o + m.X_h2s * m.m_h2s
                              + m.X_nh3 * m.m_nh3 + m.X_nh4sh * m.m_nh4sh;

        const double L_nh3    = m.r_mix * m.cp_mix * D_nh3   / m.k_mix;
        const double L_h2s    = m.r_mix * m.cp_mix * D_h2s   / m.k_mix;
        const double L_nh4sh  = m.r_mix * m.cp_mix * D_nh4sh / m.k_mix;
        const double LT_nh3   = m.r_mix * m.cp_mix * DT_nh3   / m.k_mix;
        const double LT_h2s   = m.r_mix * m.cp_mix * DT_h2s   / m.k_mix;
        const double LT_nh4sh = m.r_mix * m.cp_mix * DT_nh4sh / m.k_mix;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (int k = 1; k < m.km-1; k++) {
            for (int j = 1; j < m.jm-1; j++) {
                const double sinthe = std::max(sin(m.the.z[j]), 0.4);  // matches sinthe_min in RungeKutta_Jup.cpp
                for (int i = 1; i < m.im-1; i++) {
                    const double rm       = m.rad.z[i];
                    const double rmsinthe = rm * sinthe;

                    double dtdr = 0.0, dtdthe = 0.0, dtdphi = 0.0;
                    double dnh3dr  = 0.0, dnh3dthe  = 0.0, dnh3dphi  = 0.0;
                    double dh2sdr  = 0.0, dh2sdthe  = 0.0, dh2sdphi  = 0.0;
                    double dnh4shdr= 0.0, dnh4shdthe= 0.0, dnh4shdphi= 0.0;

                    derivative_1_order(i, j, k, dtdr,     dtdthe,     dtdphi,     m.t);
                    derivative_1_order(i, j, k, dh2sdr,   dh2sdthe,   dh2sdphi,   m.h2s);
                    derivative_1_order(i, j, k, dnh3dr,   dnh3dthe,   dnh3dphi,   m.nh3);
                    derivative_1_order(i, j, k, dnh4shdr, dnh4shdthe, dnh4shdphi, m.nh4sh);


                    derivative_1_order_boundary(i, j, k, dtdr,     dtdthe,     dtdphi,     m.t);
                    derivative_1_order_boundary(i, j, k, dh2sdr,   dh2sdthe,   dh2sdphi,   m.h2s);
                    derivative_1_order_boundary(i, j, k, dnh3dr,   dnh3dthe,   dnh3dphi,   m.nh3);
                    derivative_1_order_boundary(i, j, k, dnh4shdr, dnh4shdthe, dnh4shdphi, m.nh4sh);

                    const double dt_div = dtdr + dtdthe / rm + dtdphi / rmsinthe;
                    m.jT_nh3.x[i][j][k]   = m.r_mix * LT_nh3   / m.t.x[i][j][k] * dt_div;
                    m.jT_h2s.x[i][j][k]   = m.r_mix * LT_h2s   / m.t.x[i][j][k] * dt_div;
                    m.jT_nh4sh.x[i][j][k] = m.r_mix * LT_nh4sh / m.t.x[i][j][k] * dt_div;

                    const double dnh3_div   = dnh3dr   + dnh3dthe   / rm + dnh3dphi   / rmsinthe;
                    const double dh2s_div   = dh2sdr   + dh2sdthe   / rm + dh2sdphi   / rmsinthe;
                    const double dnh4sh_div = dnh4shdr + dnh4shdthe / rm + dnh4shdphi / rmsinthe;

                    const double cM_nh3   = L_nh3   * m.nh3.x[i][j][k]   / m.r_mix;
                    const double cM_h2s   = L_h2s   * m.h2s.x[i][j][k]   / m.r_mix;
                    const double cM_nh4sh = L_nh4sh * m.nh4sh.x[i][j][k] / m.r_mix;

                    m.j_nh3.x[i][j][k] =
                        m.m_nh3/M_mix     * (L_nh3 * dnh3_div
                        - M_mix/m.m_nh3   * cM_nh3 * dnh3_div
                        - M_mix/m.m_h2s   * cM_nh3 * dh2s_div
                        - M_mix/m.m_nh4sh * cM_nh3 * dnh4sh_div)
                        - m.jT_nh3.x[i][j][k];

                    m.j_h2s.x[i][j][k] =
                        m.m_h2s/M_mix     * (L_h2s * dh2s_div
                        - M_mix/m.m_nh3   * cM_h2s * dnh3_div
                        - M_mix/m.m_h2s   * cM_h2s * dh2s_div
                        - M_mix/m.m_nh4sh * cM_h2s * dnh4sh_div)
                        - m.jT_h2s.x[i][j][k];

                    m.j_nh4sh.x[i][j][k] =
                        m.m_nh4sh/M_mix   * (L_nh4sh * dnh4sh_div
                        - M_mix/m.m_nh3   * cM_nh4sh * dnh3_div
                        - M_mix/m.m_h2s   * cM_nh4sh * dh2s_div
                        - M_mix/m.m_nh4sh * cM_nh4sh * dnh4sh_div)
                        - m.jT_nh4sh.x[i][j][k];

                    m.difflux_h2s.x[i][j][k]   = D_h2s   * laplacian_spherical(i, j, k, m.h2s);
                    m.difflux_nh3.x[i][j][k]   = D_nh3   * laplacian_spherical(i, j, k, m.nh3);
                    m.difflux_nh4sh.x[i][j][k] = D_nh4sh * laplacian_spherical(i, j, k, m.nh4sh);

                    m.thermalmassflux.x[i][j][k] =
                          (m.j_nh3.x[i][j][k]   * m.cp_nh3
                         + m.j_h2s.x[i][j][k]   * m.cp_h2s
                         + m.j_nh4sh.x[i][j][k] * m.cp_nh4sh) * dt_div
                         + m.t.x[i][j][k]       * (m.w_nh3.x[i][j][k]   * m.k_nh3
                                                 + m.w_h2s.x[i][j][k]   * m.k_h2s
                                                 + m.w_nh4sh.x[i][j][k] * m.k_nh4sh);
                }
            }
        }

        auto end     = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for DiffMassFluxJup\n", elapsed.count() * 1e-9);
        cout << "      ATJUP: DiffMassFluxJup ended" << endl;
    }

    void ThermalPropertiesJup()
    {
        using namespace std;
        cout << endl << "      ATJUP: ThermalPropertiesJup" << endl;

        auto begin = chrono::high_resolution_clock::now();

        const double r_mixture = 1.326;                                 // [kg/m³]

        m.rg_mix = m.rg_h2 + m.rg_he + m.rg_h2s + m.rg_nh3 + m.rg_nh4sh + m.rg_h2o;
        m.r_mix  = m.r_h2  + m.r_he  + m.r_h2s  + m.r_nh3  + m.r_nh4sh  + m.r_h2o;
        m.c_mix  = m.c_h2  + m.c_he  + m.c_h2s  + m.c_nh3  + m.c_nh4sh  + m.c_h2o;

        const double M_mix = m.r_mix / m.c_mix;

        m.cp_mix  = (m.r_h2 * m.cp_h2 + m.r_h2s * m.cp_h2s
                  + m.r_he * m.cp_he  + m.r_nh3 * m.cp_nh3
                  + m.r_nh4sh * m.cp_nh4sh + m.r_h2o * m.cp_h2o) / m.r_mix;

        m.mue_mix = (m.r_h2 * m.mue_h2 + m.r_h2s * m.mue_h2s
                  + m.r_he * m.mue_he  + m.r_nh3 * m.mue_nh3
                  + m.r_nh4sh * m.mue_nh4sh + m.r_h2o * m.mue_h2o) / m.r_mix;

        m.k_mix   = (m.r_h2 * m.k_h2 + m.r_h2s * m.k_h2s
                  + m.r_he * m.k_he  + m.r_nh3 * m.k_nh3
                  + m.r_nh4sh * m.k_nh4sh + m.r_h2o * m.k_h2o) / m.r_mix;

        m.R_mix   = (m.r_h2 * m.R_h2 + m.r_he * m.R_he + m.r_h2o * m.R_h2o
                  +  m.r_h2s * m.R_h2s + m.r_nh3 * m.R_nh3 + m.r_nh4sh * m.R_nh4sh)
                  / (m.r_h2 + m.r_he + m.r_h2o + m.r_h2s + m.r_nh3 + m.r_nh4sh);

        cout.precision(10);
        cout.setf(ios::fixed);
        cout << endl
            << "     r_mixture[kg/m³] = " << r_mixture << endl << endl

            << "     rg_h2[kg/m³] = "    << m.rg_h2    << endl
            << "     rg_he[kg/m³] = "    << m.rg_he    << endl
            << "     rg_nh3[kg/m³] = "   << m.rg_nh3   << endl
            << "     rg_h2s[kg/m³] = "   << m.rg_h2s   << endl
            << "     rg_nh4sh[kg/m³] = " << m.rg_nh4sh << endl
            << "     rg_h2o[kg/m³] = "   << m.rg_h2o   << endl
            << "     rg_mix[kg/m³] = "   << m.rg_mix   << endl << endl

            << "     r_h2[kg/m³] = "    << m.r_h2    << endl
            << "     r_he[kg/m³] = "    << m.r_he    << endl
            << "     r_nh3[kg/m³] = "   << m.r_nh3   << endl
            << "     r_h2s[kg/m³] = "   << m.r_h2s   << endl
            << "     r_nh4sh[kg/m³] = " << m.r_nh4sh << endl
            << "     r_h2o[kg/m³] = "   << m.r_h2o   << endl
            << "     r_mix[kg/m³] = "   << m.r_mix   << endl << endl

            << "     c_h2[kmol/m³] = "    << m.c_h2    << endl
            << "     c_he[kmol/m³] = "    << m.c_he    << endl
            << "     c_nh3[kmol/m³] = "   << m.c_nh3   << endl
            << "     c_h2s[kmol/m³] = "   << m.c_h2s   << endl
            << "     c_nh4sh[kmol/m³] = " << m.c_nh4sh << endl
            << "     c_h2o[kmol/m³] = "   << m.c_h2o   << endl
            << "     c_mix[kmol/m³] = "   << m.c_mix   << endl << endl

            << "     ep_h2[kg/m³] = "    << m.ep_h2    << endl
            << "     ep_he[kg/m³] = "    << m.ep_he    << endl
            << "     ep_nh3[kg/m³] = "   << m.ep_nh3   << endl
            << "     ep_h2s[kg/m³] = "   << m.ep_h2s   << endl
            << "     ep_nh4sh[kg/m³] = " << m.ep_nh4sh << endl
            << "     ep_h2o[kg/m³] = "   << m.ep_h2o   << endl << endl

            << "     X_h2[/] = "    << m.X_h2    << endl
            << "     X_he[/] = "    << m.X_he    << endl
            << "     X_nh3[/] = "   << m.X_nh3   << endl
            << "     X_h2s[/] = "   << m.X_h2s   << endl
            << "     X_nh4sh[/] = " << m.X_nh4sh << endl
            << "     X_h2o[/] = "   << m.X_h2o   << endl << endl

            << "     gam[/] = "        << m.gam    << endl
            << "     M_mix[kg/Kmol] = " << M_mix   << endl
            << "     r_mix[kg/m³] = "  << m.r_mix  << endl
            << "     c_mix[kmol/m³] = " << m.c_mix << endl << endl

            << "     cp_mix[J/(kg*K)] = "  << m.cp_mix  << endl
            << "     mue_mix[Ns/m²] = "    << m.mue_mix << endl
            << "     k_mix[W/(m*K)] = "    << m.k_mix   << endl
            << "     R_mix[J/(kg*K)] = "   << m.R_mix   << endl << endl;

        auto end     = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ThermalPropertiesJup\n", elapsed.count() * 1e-9);
        cout << "      ATJUP: ThermalPropertiesJup ended" << endl;
    }

private:
    cJupiterModel& m;

    static double react_rate_const(double T_K, double T_d, double A, double B)
    {
        return A * pow(T_K, B) * exp(-T_d / T_K);
    }

    void derivative_1_order(int i, int j, int k,
        double& dcdr, double& dcdthe, double& dcdphi, Array& c)
    {
        const double dr_2   = 2.0 * m.dr;
        const double dthe_2 = 2.0 * m.dthe;
        const double dphi_2 = 2.0 * m.dphi;
        const double exp_rm = m.coord_stretching ? 1.0 / (m.rad.z[i] + 1.0) : 1.0;

        // r-direction: cubic extrapolation at solid surface, central diff at adjacent fluid
        if ((m.SeaMount.x[i][j][k] == 1.0) && (m.SeaMount.x[i+1][j][k] == 0.0)
            && (m.SeaMount.x[i+2][j][k] == 0.0)) {
            c.x[i][j][k] = c.x[i+3][j][k] - 3.0*c.x[i+2][j][k] + 3.0*c.x[i+1][j][k];
            dcdr = (-3.0*c.x[i][j][k] + 4.0*c.x[i+1][j][k] - c.x[i+2][j][k]) / dr_2;
        }
        if ((m.SeaMount.x[i-1][j][k] == 1.0) && (m.SeaMount.x[i][j][k] == 0.0)
            && (m.SeaMount.x[i+1][j][k] == 0.0))
            dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k]) / dr_2;
        if ((m.SeaMount.x[i-1][j][k] == 0.0) && (m.SeaMount.x[i][j][k] == 0.0)
            && (m.SeaMount.x[i+1][j][k] == 1.0))
            dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k]) / dr_2;
        if ((m.SeaMount.x[i-1][j][k] == 0.0) && (m.SeaMount.x[i][j][k] == 0.0)
            && (m.SeaMount.x[i+1][j][k] == 0.0))
            dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k]) / dr_2;
        dcdr *= exp_rm;

        // theta-direction: cubic extrapolation at solid surface, central diff at adjacent fluid
        if ((m.SeaMount.x[i][j][k] == 1.0) && (m.SeaMount.x[i][j+1][k] == 0.0)
            && (m.SeaMount.x[i][j+2][k] == 0.0)) {
            c.x[i][j][k] = 4.0/3.0*c.x[i][j+1][k] - 1.0/3.0*c.x[i][j+2][k];
            dcdthe = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j+1][k] - c.x[i][j+2][k]) / dthe_2;
        }
        if ((m.SeaMount.x[i][j-1][k] == 1.0) && (m.SeaMount.x[i][j][k] == 0.0)
            && (m.SeaMount.x[i][j+1][k] == 0.0))
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / dthe_2;
        if ((m.SeaMount.x[i][j-1][k] == 0.0) && (m.SeaMount.x[i][j][k] == 0.0)
            && (m.SeaMount.x[i][j+1][k] == 1.0))
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / dthe_2;

        if ((m.SeaMount.x[i][j][k] == 1.0) && (m.SeaMount.x[i][j-1][k] == 0.0)
            && (m.SeaMount.x[i][j-2][k] == 0.0)) {
            c.x[i][j][k] = 4.0/3.0*c.x[i][j-1][k] - 1.0/3.0*c.x[i][j-2][k];
            dcdthe = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j-1][k] - c.x[i][j-2][k]) / dthe_2;
        }
        if ((m.SeaMount.x[i][j-1][k] == 0.0) && (m.SeaMount.x[i][j][k] == 0.0)
            && (m.SeaMount.x[i][j+1][k] == 0.0))
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / dthe_2;

        // phi-direction: cubic extrapolation at solid surface, central diff at adjacent fluid
        if ((m.SeaMount.x[i][j][k] == 1.0) && (m.SeaMount.x[i][j][k+1] == 0.0)
            && (m.SeaMount.x[i][j][k+2] == 0.0)) {
            c.x[i][j][k] = 4.0/3.0*c.x[i][j][k+1] - 1.0/3.0*c.x[i][j][k+2];
            dcdphi = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j][k+1] - c.x[i][j][k+2]) / dphi_2;
        }
        if ((m.SeaMount.x[i][j][k-1] == 1.0) && (m.SeaMount.x[i][j][k] == 0.0)
            && (m.SeaMount.x[i][j][k+1] == 0.0))
            dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1]) / dphi_2;
        if ((m.SeaMount.x[i][j][k-1] == 0.0) && (m.SeaMount.x[i][j][k] == 0.0)
            && (m.SeaMount.x[i][j][k+1] == 1.0))
            dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1]) / dphi_2;

        if ((m.SeaMount.x[i][j][k] == 1.0) && (m.SeaMount.x[i][j][k-1] == 0.0)
            && (m.SeaMount.x[i][j][k-2] == 0.0)) {
            c.x[i][j][k] = 4.0/3.0*c.x[i][j][k-1] - 1.0/3.0*c.x[i][j][k-2];
            dcdphi = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j][k-1] - c.x[i][j][k-2]) / dphi_2;
        }
        if ((m.SeaMount.x[i][j][k-1] == 0.0) && (m.SeaMount.x[i][j][k] == 0.0)
            && (m.SeaMount.x[i][j][k+1] == 0.0))
            dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1]) / dphi_2;
    }

    void derivative_1_order_boundary(int i, int j, int k,
        double& dcdr, double& dcdthe, double& dcdphi, Array& c)
    {
        if (i == 0) {
            const double exp_rm = m.coord_stretching ? 1.0 / (m.rad.z[0] + 1.0) : 1.0;
            dcdr = (-3.0*c.x[0][j][k]      + 4.0*c.x[1][j][k]      - c.x[2][j][k])      / (2.0*m.dr) * exp_rm;
        }
        if (i == m.im-1) {
            const double exp_rm = m.coord_stretching ? 1.0 / (m.rad.z[m.im-1] + 1.0) : 1.0;
            dcdr = (-3.0*c.x[m.im-1][j][k] + 4.0*c.x[m.im-2][j][k] - c.x[m.im-3][j][k]) / (2.0*m.dr) * exp_rm;
        }

        if (j == 0)
            dcdthe = (-3.0*c.x[i][0][k]       + 4.0*c.x[i][1][k]        - c.x[i][2][k])        / (2.0*m.dthe);
        if (j == m.jm-1)
            dcdthe = (-3.0*c.x[i][m.jm-1][k]  + 4.0*c.x[i][m.jm-2][k]  - c.x[i][m.jm-3][k])  / (2.0*m.dthe);

        if (k == 0)
            dcdphi = (-3.0*c.x[i][j][0]       + 4.0*c.x[i][j][1]        - c.x[i][j][2])        / (2.0*m.dphi);
        if (k == m.km-1)
            dcdphi = (-3.0*c.x[i][j][m.km-1]  + 4.0*c.x[i][j][m.km-2]  - c.x[i][j][m.km-3])  / (2.0*m.dphi);
    }

    // Full spherical Laplacian: ∂²c/∂r² + (2/r)∂c/∂r
    //   + (1/r²)[∂²c/∂θ² + (cosθ/sinθ)∂c/∂θ]
    //   + (1/(r sinθ))² ∂²c/∂φ²
    // exp_rm/exp_2_rm scale model-coordinate r-derivatives to physical ones.
    double laplacian_spherical(int i, int j, int k, Array& c)
    {
        const double rm       = m.rad.z[i];
        const double exp_rm   = m.coord_stretching ? 1.0 / (rm + 1.0) : 1.0;
        const double exp_2_rm = exp_rm * exp_rm;
        const double sinthe   = sin(m.the.z[j]);
        const double costhe   = cos(m.the.z[j]);
        const double rmsinthe = rm * std::max(sinthe, 0.4);  // matches sinthe_min in RungeKutta_Jup.cpp

        const double d2cdr2   = (c.x[i+1][j][k] - 2.0*c.x[i][j][k] + c.x[i-1][j][k]) / (m.dr   * m.dr) * exp_2_rm;
        const double d2cdthe2 = (c.x[i][j+1][k] - 2.0*c.x[i][j][k] + c.x[i][j-1][k]) / (m.dthe * m.dthe);
        const double d2cdphi2 = (c.x[i][j][k+1] - 2.0*c.x[i][j][k] + c.x[i][j][k-1]) / (m.dphi * m.dphi);

        const double dcdr   = (c.x[i+1][j][k] - c.x[i-1][j][k]) / (2.0 * m.dr) * exp_rm;
        const double dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / (2.0 * m.dthe);

        return d2cdr2
             + 2.0 / rm * dcdr
             + d2cdthe2 / (rm * rm)
             + costhe / (rm * rmsinthe) * dcdthe
             + d2cdphi2 / (rmsinthe * rmsinthe);
    }
};
