/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * Dry convective adjustment.
 *
 * Reference algorithm:
 *   Manabe, S. and Strickler, R. F.: "Thermal Equilibrium of the Atmosphere with a Convective
 *   Adjustment", J. Atmos. Sci. 21, 361-385, 1964.
 *
 * WHY THE MODEL NEEDS THIS. Measured from its own restart files against its own cp_mix, ATJUP
 * develops a thin superadiabatic layer and deepens it as a run proceeds: at iteration 100 two
 * levels around 21-24 km exceed the dry adiabat by 0.016 K/km, at iteration 200 by 0.035, and by
 * iteration 500 five levels from 21 to 35 km exceed it by up to 0.125. It can accumulate that
 * because the momentum equation never felt the buoyancy — with the body forces in the wrong unit
 * system nothing responded to the instability, so nothing relieved it.
 *
 * Once ATJUP_NONDIM switches the buoyancy on, the vertical velocity grows steadily (max|u| from 38
 * to 323 m/s over 250 iterations, e-folding about 143), and the adiabatic cooling of the runaway
 * updraft drives the coldest cell from 110 K to 20 K and then through zero: the run goes
 * non-finite at iteration 233. How much of that growth is this instability is not settled — the
 * layer is thin and weak, and the e-folding its N implies is 207 iterations rather than the 143
 * observed, so the horizontal equator-to-pole density contrast is very likely the larger driver.
 * This class removes one of the two, which is a precondition for judging the other.
 *
 * The buoyancy term cannot fix this by itself, and it is worth being clear about why. It is written
 * as an anomaly about the horizontal mean of each level, so a purely one-dimensional superadiabatic
 * column produces exactly zero force. The term responds to horizontal density contrasts; the
 * unstable stratification is what makes those contrasts grow rather than oscillate. Nothing in the
 * model restores a column to its adiabat, which is what this class does.
 *
 * WHAT IT DOES. Each column is swept from the bottom up. Wherever a layer pair is steeper than the
 * dry adiabat, the pair is mixed to exactly the adiabatic lapse rate, conserving the mass-weighted
 * enthalpy of the pair. Sweeps repeat until the column is stable. This is the classical scheme, and
 * its two properties are the ones that matter here: it removes the instability completely rather
 * than damping it, and it does not create or destroy energy.
 *
 * WHAT IT DELIBERATELY DOES NOT DO. It mixes temperature only, not composition. Real convection
 * carries the species with it, and a moist adjustment would use the saturated adiabat where cloud
 * is present rather than the dry one. Both are reasonable extensions; neither is done here, because
 * each is a modelling decision with its own consequences for the microphysics that already runs in
 * SaturationAdjustmentJup and PrecipitationJup.
 */

#pragma once

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

class cJupiterModel;

class ConvectiveAdjustmentJup {
public:
    explicit ConvectiveAdjustmentJup(cJupiterModel& model) : m(model) {}

    void run();

private:
    cJupiterModel& m;
};
