/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * ATJUP_CWB_DIAG=1 -- the column water budget of H2O (vapour + cloud + ice), PRINT-ONLY, default off.
 *
 * WHY (2026-10-08). The saturation adjustment's own budget (ATJUP_SATADJ_DIAG) showed that the H2O
 * column loses 3.3 % in 224 iterations and that the adjustment is not where it goes: a steady
 * -680 g/m2 per pair of iterations disappears BETWEEN its calls, 30 000 times the printed surface
 * precipitation. This instrument marks the column after every stage of the time loop that can
 * write the three fields and charges the difference to that stage. Mirrored from
 * ATOM_Precipitation's ATM_CWB_DIAG.
 *
 * The Runge-Kutta step is split further, because it is two things at once:
 *   rk_reset   it integrates from the n-level copies (h2on, ...), which restoreVar refreshed at
 *              the END of the previous iteration -- so every direct write made since then (the
 *              saturation adjustment, an in-place precipitation depletion) is discarded here;
 *   rk terms   the RK4-weighted column integrals of each right-hand-side term of the three
 *              equations, accumulated per thread inside RHSJup: radial / meridional / zonal
 *              transport (advective form, -u.grad q), diffusion, the precipitation source;
 *   rk_rest    what the step changed beyond reset + terms: the max(0, .) floor of the final
 *              update and round-off. An identity check when the floor does not bind.
 * q*div(u) is accumulated beside them and NOT applied: the transport is written in advective
 * form, and -u.grad q = -div(q u) + q div(u), so on a closed column the advective form creates or
 * destroys exactly the integral of q*div(u). Printed so the two can be compared.
 * With ATJUP_SPECIES_DIVU set the term IS applied (see RHS_Jup_Turb.cpp) and appears as its own
 * row, "flux-form term -s q div(u)"; the not-applied row stays as the reference.
 *
 * Unit: g/m2, the mean over the sphere weighted by sin(colatitude), layer thickness in metres on
 * densities in kg/m3 -- the same sum the saturation adjustment's instrument prints. The top level
 * has no thickness and is not counted.
 *
 * The state lives in a function-local static so that the three translation units share one
 * instance without a model member (C++11 has no inline variables).
 */

#pragma once

namespace JupCwb {

    enum { T_RAD = 0, T_THE, T_PHI, T_DIFF, T_PRECIP, T_DIVU, T_QDIV, NTERM };
    const int MAXTHR = 256;

    struct State {
        bool   on;                       // ATJUP_CWB_DIAG
        double stage_wgt;                // RK4 weight of the stage being evaluated (1, 2, 2, 1)
        double acc[MAXTHR][8];           // [thread][term], weighted sums of the tendencies
    };

    inline State& st(){
        static State s = State();
        return s;
    }
}
