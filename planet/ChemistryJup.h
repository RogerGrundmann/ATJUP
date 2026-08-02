#pragma once

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>   // getenv/atoi for the NH4SH rate-law knobs
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
#include "ATPhys.h"
#include "FluxLimiter.h"

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

        // ====================================================================================
        // TWO DEFECTS IN THE RATE LAW BELOW, EACH WITH A KNOB THAT RESTORES THE OLD BEHAVIOUR.
        //
        // What they cost, measured: two runs of 450 iterations at dt=0.001, identical except
        // <chemical_reaction> 0 against 1. On the meridional slice k=180 the whole NH4SH field
        // is the chemistry's doing — with the switch off it is not merely smaller but exactly
        // zero everywhere, so neither transport nor Stokes settling contributes anything:
        //
        //     i=12 j=9   (81N)   245.4  ->  0        i=13 j=151 (61S)   242.7  ->  0
        //     i=12 j=8   (82N)   205.2  ->  0        i=12 j=152 (62S)   151.1  ->  0
        //
        // And the balance does not close: producing 245.4 units of NH4SH moved NH3 by 1e-4
        // (0.0067 against 0.0068) and the temperature not at all. Six orders of magnitude.
        // The vertical profile at j=9 is a geometric series, factor 4-6 per level over ten
        // levels, peaking at i=12 — which is BELOW the 200..230 K reaction layer (i=14..19
        // there), i.e. where the gate is shut.
        //
        // (1) ATJUP_CHEM_GATE_ZERO — the temperature gate had no else branch. w_nh3/w_h2s/
        //     w_nh4sh are only assigned inside `if (t_00 <= t_u <= t_0)`; they are persistent
        //     Arrays, so a cell that LEAVES the window keeps the last rate it ever had and goes
        //     on applying it to the tendency every iteration for the rest of the run. All the
        //     hot spots sit outside the window: (12,9) 243.5 K, (12,151) 240.0 K, (13,151)
        //     236.2 K. Default 1 zeroes the three rates outside the window; 0 is the old freeze.
        //
        // (2) ATJUP_CHEM_MOLAR_CONC — DEFAULT 0, and the default is a modelling decision, not an
        //     endorsement of the old formula. Measured, 450 iterations at dt=0.001, otherwise
        //     identical runs, NH4SH on the meridional slice at iteration 450:
        //
        //         old normalisation   i=12 j=9  245.5    i=13 j=151  242.7    max 245.5
        //         (1) alone            "        12.38     "           33.49   max  82.08
        //         (1) + this           "         0        "            0      max   0
        //
        //     The repair is dimensionally right and it removes the NH4SH cloud completely (below
        //     the 1e-6 output precision). The arithmetic says why: the old denominator inflated
        //     r_mix = 0.42 to r_mix/sum_c ~ 1000, a factor ~2400 per concentration and ~6e6 in
        //     the quadratic kf*c_nh3*c_h2s, so switching to the true concentration collapses the
        //     forward rate by those six orders. A = 15000 and T_d = 3020 above were evidently
        //     calibrated against the inflated values; correcting the concentrations without
        //     recalibrating the rate leaves no ammonium hydrosulfide at all. Until the rate is
        //     re-fitted, the shipped default keeps the old normalisation and relies on (1),
        //     which lowers the hot spots by a factor 20 and touches no calibration.
        //
        //     What it does when set to 1: the concentrations were renormalised by
        //         denom = r_mix / sum_c,   sum_c = nh3/m_nh3 + h2s/m_h2s + nh4sh/m_nh4sh
        //     which forces c_nh3 + c_h2s + c_nh4sh = r_mix: three trace species rescaled to
        //     carry the ENTIRE mixture density. Measured at the hot spots, NH4SH then holds
        //     99.99 % of sum_c (75 % already at NH4SH = 0.064), so c_nh4sh = r_mix to four
        //     digits no matter how large the field grows — the back reaction kb*c_nh4sh
        //     saturates and the rate law stops seeing its own product. It is also a unit
        //     mismatch: nh3 and h2s are mass fractions [kg/kg] while the NH4SH initial field is
        //     built as CORR_NH4SH * r_mix * q_Rain, a density [kg/m3], and sum_c adds them.
        //     Set to 1 it uses the plain molar concentration c_x = rho_mix * w_x / m_x
        //     [kmol/m3]; 0 (the default, see above) keeps the normalisation.
        //
        // rho_mix is refreshed by computeMixtureDensity() AFTER this routine in the physics
        // block (cJupiterModel.cpp), so on the first iteration it is still empty; fall back to
        // the scalar r_mix then, as the RHS does for the same reason.
        // ====================================================================================
        static const int gate_zero  = [](){
            const char* e = getenv("ATJUP_CHEM_GATE_ZERO");  return e ? atoi(e) : 1; }();
        static const int molar_conc = [](){
            const char* e = getenv("ATJUP_CHEM_MOLAR_CONC"); return e ? atoi(e) : 0; }();

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
                        double denom;
                        if (molar_conc) {
                            // c_x = rho * w_x / m_x, the plain molar concentration [kmol/m3].
                            double rho = m.rho_mix.x[i][j][k];
                            if (!(rho > 0.0)) rho = m.r_mix;   // first iteration: not filled yet
                            denom = rho;
                        } else {
                            const double sum_c = m.nh3.x[i][j][k]   / m.m_nh3
                                               + m.h2s.x[i][j][k]   / m.m_h2s
                                               + m.nh4sh.x[i][j][k] / m.m_nh4sh;
                            denom = (sum_c == 0.0) ? 0.0 : m.r_mix / sum_c;
                        }

                        const double c_nh3   = m.nh3.x[i][j][k]   / m.m_nh3   * denom;
                        const double c_h2s   = m.h2s.x[i][j][k]   / m.m_h2s   * denom;
                        const double c_nh4sh = m.nh4sh.x[i][j][k] / m.m_nh4sh * denom;

                        const double R_diff = kf * c_nh3 * c_h2s - kb * c_nh4sh;

                        m.w_nh3.x[i][j][k]   = -m.m_nh3   * R_diff;
                        m.w_h2s.x[i][j][k]   = -m.m_h2s   * R_diff;
                        m.w_nh4sh.x[i][j][k] =  m.m_nh4sh * R_diff;
                    } else if (gate_zero) {
                        // Outside the 200..230 K formation window there is no reaction, so the
                        // rates are zero — not "whatever they were the last time this cell was
                        // inside the window".
                        m.w_nh3.x[i][j][k]   = 0.0;
                        m.w_h2s.x[i][j][k]   = 0.0;
                        m.w_nh4sh.x[i][j][k] = 0.0;
                    }

                    // PLUS, not minus. massflux_* is added to the species tendency in
                    // RHS_Jup_Turb.cpp, and difflux_* is computed just above as
                    // D_x * laplacian(c_x) — the diffusive TENDENCY, positive where a species
                    // sits in a local minimum, which is exactly what Fick gives:
                    //     dc/dt = -div(j) + sources,   j = -D grad(c),   -div(j) = +D lap(c)
                    // Subtracting it made the multicomponent diffusion an ANTI-diffusion,
                    // sharpening every species gradient instead of smoothing it.
                    //
                    // It survived because the coefficient is tiny: D_x is built as
                    // mue_x/(rg_x*sc_x) with rg_h2s = 949 kg/m3, the density of the CONDENSED
                    // phase rather than the gas, giving D = 1.4e-8 m2/s, and diff_* prints as
                    // 0.000000 kg/(m3 s) against massflux_* of 3e-4.
                    //
                    // Tiny is NOT the same as harmless, and the measurement says so. Anti-
                    // diffusion is self-amplifying: it sharpens a gradient, the sharper gradient
                    // raises the Laplacian, which sharpens it further. Differencing the full 3D
                    // state at iteration 100 against the uncorrected run, everything else equal:
                    //
                    //   nh3    1.27 %      nh4sh  1.56 %      h2s   0.24 %
                    //   t      0.027 %     w      0.066 %
                    //
                    // So a term that never shows up in printMinMax had moved the ammonia field
                    // by more than a percent in 100 iterations. The questionable D_x is left
                    // alone here — it belongs to a separate question.
                    //
                    // Also left alone, but worth recording: `chemical_reaction` in the RHS
                    // multiplies the WHOLE of massflux_*, so setting that switch to 0 to disable
                    // the chemistry silently disables this diffusion as well.
                    m.massflux_h2s.x[i][j][k]   = m.w_h2s.x[i][j][k]   + m.difflux_h2s.x[i][j][k];
                    m.massflux_nh3.x[i][j][k]   = m.w_nh3.x[i][j][k]   + m.difflux_nh3.x[i][j][k];
                    m.massflux_nh4sh.x[i][j][k] = m.w_nh4sh.x[i][j][k] + m.difflux_nh4sh.x[i][j][k];
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

//        const double h2   = 0.74;
//        const double he   = 0.24;

        const double DT_nh3   = m.mue_nh3   / (m.rg_nh3   * m.sc_nh3);
        const double DT_h2s   = m.mue_h2s   / (m.rg_h2s   * m.sc_h2s);
        const double DT_nh4sh = m.mue_nh4sh / (m.rg_nh4sh * m.sc_nh4sh);
        const double D_nh3    = DT_nh3;
        const double D_h2s    = DT_h2s;
        const double D_nh4sh  = DT_nh4sh;

        const double M_mix    = m.X_h2 * m.m_h2 + m.X_he * m.m_he
                              + m.X_h2o * m.m_h2o + m.X_h2s * m.m_h2s
                              + m.X_nh3 * m.m_nh3 + m.X_nh4sh * m.m_nh4sh;

        const double L_nh3    = m.r_mix * m.cp_mix * D_nh3    / m.k_mix;
        const double L_h2s    = m.r_mix * m.cp_mix * D_h2s    / m.k_mix;
        const double L_nh4sh  = m.r_mix * m.cp_mix * D_nh4sh  / m.k_mix;
        const double LT_nh3   = m.r_mix * m.cp_mix * DT_nh3   / m.k_mix;
        const double LT_h2s   = m.r_mix * m.cp_mix * DT_h2s   / m.k_mix;
        const double LT_nh4sh = m.r_mix * m.cp_mix * DT_nh4sh / m.k_mix;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (int k = 1; k < m.km-1; k++) {
            for (int j = 1; j < m.jm-1; j++) {
                const double sinthe = std::max(sin(m.the.z[j]),
                                               ATPhys::polar_divisor_floor<cJupiterModel>());
                for (int i = 1; i < m.im-1; i++) {
                    const double rm       = m.rad.z[i];
                    const double rmsinthe = rm * sinthe;

                    // The thermo-diffusion fluxes below divide by the local temperature, which
                    // bcSolidGround sets to exactly 0 inside the SeaMount — 1/0 gave inf (or
                    // NaN once dt_div was 0 as well) in every solid cell, and the inf then
                    // spread through massflux_* into the species equations. Solid cells carry
                    // no diffusive flux, so zero them and move on. !(t > 0.0) rather than
                    // (t <= 0.0) so an incoming NaN is caught rather than propagated.
                    if(!(m.t.x[i][j][k] > 0.0) || m.SeaMount.x[i][j][k] == 1.0){
                        m.jT_nh3.x[i][j][k]   = 0.0;
                        m.jT_h2s.x[i][j][k]   = 0.0;
                        m.jT_nh4sh.x[i][j][k] = 0.0;
                        m.j_nh3.x[i][j][k]    = 0.0;
                        m.j_h2s.x[i][j][k]    = 0.0;
                        m.j_nh4sh.x[i][j][k]  = 0.0;
                        m.difflux_h2s.x[i][j][k]   = 0.0;
                        m.difflux_nh3.x[i][j][k]   = 0.0;
                        m.difflux_nh4sh.x[i][j][k] = 0.0;
                        m.thermalmassflux.x[i][j][k] = 0.0;
                        continue;
                    }

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

                    // Pole-symmetric sum of gradient components: dXdr and dXdphi are
                    // pole-symmetric for symmetric inputs, dXdthe is pole-antisymmetric.
                    // std::abs() on the theta term symmetrizes the resulting scalar.
                    // (difflux_* below uses laplacian_spherical() which is already
                    // pole-symmetric by construction via cosθ·∂c/∂θ; no fix needed there.)
                    const double dt_div = dtdr + std::abs(dtdthe) / rm + dtdphi / rmsinthe;

                    // ATJUP_LOCAL_RHO: the Lewis groups L = rho*cp*D/k and LT = rho*cp*DT/k are
                    // precomputed once with the constant r_mix, so rescaling them per cell is a
                    // single factor rho_c/r_mix. Where that factor lands is worth spelling out,
                    // because it is not uniform:
                    //   jT_*  carries rho TWICE (once explicitly, once inside LT) -> (rho/r_mix)^2
                    //   cM_*  divides by rho again, so the two cancel and it is UNCHANGED
                    //   j_*   the bare L terms scale linearly with rho/r_mix
                    // With rho_c/r_mix reaching 1/200 near the top, the thermo-diffusion flux
                    // there is cut by ~4e-5 while the concentration part is only cut by 1/200.
                    const double rho_c = m.rho_at(i, j, k);
                    const double sc    = rho_c / m.r_mix;      // 1.0 exactly when the knob is off

                    // ONE density, not two. This line used to read
                    //     jT = m.r_mix * LT_x / t * dt_div
                    // with LT_x = r_mix * cp_mix * DT_x / k_mix, so the mixture density entered
                    // TWICE and the flux went as rho^2. The Soret (thermal diffusion) mass flux
                    // is j_T = -rho * D_T * grad(T)/T — linear in the density, once. The double
                    // count was invisible while the density was the constant 1.2844: it merely
                    // scaled every jT_* by that factor. With the local density it is not
                    // invisible at all, because rho spans 0.006 to 1.09 across the shell and the
                    // square would cut the flux near the top by 4e-5 instead of 1/200.
                    // LT_x * sc already carries exactly one density — the local one when
                    // ATJUP_LOCAL_RHO is set, r_mix otherwise — so no explicit factor here.
                    m.jT_nh3.x[i][j][k]   = (LT_nh3   * sc) / m.t.x[i][j][k] * dt_div;
                    m.jT_h2s.x[i][j][k]   = (LT_h2s   * sc) / m.t.x[i][j][k] * dt_div;
                    m.jT_nh4sh.x[i][j][k] = (LT_nh4sh * sc) / m.t.x[i][j][k] * dt_div;

                    const double dnh3_div   = dnh3dr   + std::abs(dnh3dthe)   / rm + dnh3dphi   / rmsinthe;
                    const double dh2s_div   = dh2sdr   + std::abs(dh2sdthe)   / rm + dh2sdphi   / rmsinthe;
                    const double dnh4sh_div = dnh4shdr + std::abs(dnh4shdthe) / rm + dnh4shdphi / rmsinthe;

                    const double L_nh3_c   = L_nh3   * sc;
                    const double L_h2s_c   = L_h2s   * sc;
                    const double L_nh4sh_c = L_nh4sh * sc;

                    const double cM_nh3   = L_nh3_c   * m.nh3.x[i][j][k]   / rho_c;
                    const double cM_h2s   = L_h2s_c   * m.h2s.x[i][j][k]   / rho_c;
                    const double cM_nh4sh = L_nh4sh_c * m.nh4sh.x[i][j][k] / rho_c;

                    m.j_nh3.x[i][j][k] =
                        m.m_nh3/M_mix     * (L_nh3_c * dnh3_div
                        - M_mix/m.m_nh3   * cM_nh3 * dnh3_div
                        - M_mix/m.m_h2s   * cM_nh3 * dh2s_div
                        - M_mix/m.m_nh4sh * cM_nh3 * dnh4sh_div)
                        - m.jT_nh3.x[i][j][k];

                    m.j_h2s.x[i][j][k] =
                        m.m_h2s/M_mix     * (L_h2s_c * dh2s_div
                        - M_mix/m.m_nh3   * cM_h2s * dnh3_div
                        - M_mix/m.m_h2s   * cM_h2s * dh2s_div
                        - M_mix/m.m_nh4sh * cM_h2s * dnh4sh_div)
                        - m.jT_h2s.x[i][j][k];

                    m.j_nh4sh.x[i][j][k] =
                        m.m_nh4sh/M_mix   * (L_nh4sh_c * dnh4sh_div
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

    // The TVD limiter is the SHARED FluxLimiter<Planet> template — 97 % identical between
    // the two models, the whole difference being metricRadius(). See FluxLimiter.h for why
    // this one routine could be shared while the rest of the chemistry cannot.
    void FluxLimiterNH4SH(){ FluxLimiter<cJupiterModel>(m).nh4sh(); }

    void ThermalPropertiesJup()
    {
        using namespace std;
        cout << endl << "      ATJUP: ThermalPropertiesJup" << endl;

        auto begin = chrono::high_resolution_clock::now();

        const double r_mixture = 1.326;                                 // [kg/m³]

        m.rg_mix = m.rg_h2 + m.rg_he + m.rg_h2s + m.rg_nh3 + m.rg_nh4sh + m.rg_h2o + m.rg_ch4;
        m.r_mix  = m.r_h2  + m.r_he  + m.r_h2s  + m.r_nh3  + m.r_nh4sh  + m.r_h2o  + m.r_ch4;
        m.c_mix  = m.c_h2  + m.c_he  + m.c_h2s  + m.c_nh3  + m.c_nh4sh  + m.c_h2o  + m.c_ch4;

        const double M_mix = m.r_mix / m.c_mix;

        m.cp_mix  = (m.r_h2 * m.cp_h2 + m.r_h2s * m.cp_h2s
                  + m.r_he * m.cp_he  + m.r_nh3 * m.cp_nh3
                  + m.r_nh4sh * m.cp_nh4sh + m.r_h2o * m.cp_h2o
                  + m.r_ch4 * m.cp_ch4) / m.r_mix;

        m.mue_mix = (m.r_h2 * m.mue_h2 + m.r_h2s * m.mue_h2s
                  + m.r_he * m.mue_he  + m.r_nh3 * m.mue_nh3
                  + m.r_nh4sh * m.mue_nh4sh + m.r_h2o * m.mue_h2o
                  + m.r_ch4 * m.mue_ch4) / m.r_mix;

        m.k_mix   = (m.r_h2 * m.k_h2 + m.r_h2s * m.k_h2s
                  + m.r_he * m.k_he  + m.r_nh3 * m.k_nh3
                  + m.r_nh4sh * m.k_nh4sh + m.r_h2o * m.k_h2o
                  + m.r_ch4 * m.k_ch4) / m.r_mix;

        m.R_mix   = (m.r_h2 * m.R_h2 + m.r_he * m.R_he + m.r_h2o * m.R_h2o
                  +  m.r_h2s * m.R_h2s + m.r_nh3 * m.R_nh3 + m.r_nh4sh * m.R_nh4sh
                  +  m.r_ch4 * m.R_ch4)
                  / (m.r_h2 + m.r_he + m.r_h2o + m.r_h2s + m.r_nh3 + m.r_nh4sh + m.r_ch4);
/*
        R_mix[i][j][k] = h2.x[i][j][k]  * R_h2                            // AI
                       + he.x[i][j][k]  * R_he 
                       + nh3.x[i][j][k] * R_nh3 
                       + h2s.x[i][j][k] * R_h2s 
                       + h2o.x[i][j][k] * R_h2o;
*/
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
        const double rmsinthe = rm * std::max(sinthe,
                                    ATPhys::polar_divisor_floor<cJupiterModel>());

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

    // Superbee: most compressive TVD limiter — best for sharp cloud fronts.

    // Van Leer: smooth, differentiable — good general-purpose alternative.
};
