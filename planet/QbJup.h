/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * ATJUP_QB_DIAG=<n> -- the budget of one species group by level, PRINT-ONLY, default off; printed
 * every n iterations. ATJUP_QB_SPECIES picks the group: nh3 (default; vapour + cloud + ice), h2s,
 * nh4sh.
 *
 * WHY (2026-10-19). In the default runs the species at i = 28 grow away from the obstacle (GRS_1000,
 * 22 S: NH3 ice 0.46 -> 0.96, NH3 vapour 0.0085 -> 0.07, H2S +50 % over 900 iterations) while the
 * columns standing on the obstacle top keep their start values. ATJUP_CWB_DIAG sums the H2O column
 * and cannot say which level gains at whose expense; this instrument does, for the other species.
 *
 * Same scheme as ATJUP_TB_DIAG (TbJup.h): the level mean (weighted by sin(colatitude), in g/m3, over
 * ALL cells of a level, the solid ones of the SeaMount included) is marked after every stage of the
 * time loop and the difference charged to that stage; the Runge-Kutta step is split into
 *   rk_reset   the step integrates from the n-level copies, so every direct write made since the end
 *              of the previous iteration and not copied into them is discarded here;
 *   rk terms   the RK4-weighted level means of each term of the group's right-hand sides,
 *              accumulated per thread in RHSJup: the three transport terms, diffusion split into its
 *              molecular (1/(sc re)) and turbulent (nue_t) parts, the chemistry (NH4SH formation),
 *              the precipitation source, the NH4SH settling (Stokes + flux limiter);
 *   rk_rest    what the step changed beyond reset + terms: the floor of the update and the rows and
 *              planes the step does not integrate.
 *
 * The state lives in a function-local static so that the three translation units share one
 * instance without a model member (C++11 has no inline variables).
 */

#pragma once

#include <cstdlib>
#include <cstring>

namespace JupQb {

    enum { T_TR_RAD = 0, T_TR_THE, T_TR_PHI, T_DIFF_MOL, T_DIFF_TURB, T_CHEM, T_PRECIP, T_SETTLE, NTERM };
    enum { G_NH3 = 0, G_H2S, G_NH4SH };
    const int MAXTHR = 256;
    const int MAXLEV = 96;

    struct State {
        bool   on;                               // ATJUP_QB_DIAG
        int    group;                            // ATJUP_QB_SPECIES
        double stage_wgt;                        // RK4 weight of the stage being evaluated (1, 2, 2, 1)
        double acc[MAXTHR][NTERM][MAXLEV];       // [thread][term][level], weighted sums of the tendencies
    };

    inline State& st(){
        static State s = State();
        return s;
    }

    inline int group_from_env(){
        const char* e = getenv("ATJUP_QB_SPECIES");
        if(!e) return G_NH3;
        if(!strcmp(e, "h2s"))   return G_H2S;
        if(!strcmp(e, "nh4sh")) return G_NH4SH;
        return G_NH3;
    }
}
