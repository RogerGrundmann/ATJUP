/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * Thermal-wind initialisation of the horizontal velocity field.
 *
 * WHY. With the body forces in their proper units (ATJUP_NONDIM) the model accelerates from 38 to
 * 323 m/s in 250 iterations and dies. The driver is not the vertical stratification — a convective
 * adjustment removes that and changes nothing — but the HORIZONTAL density contrast: measured per
 * level, the spread of rho about its own level mean implies a buoyancy of 8 to 10 in the units
 * rhs_u works in, against transport terms of order one.
 *
 * That contrast is not itself the problem. A buoyancy field B(r,theta,phi) in the radial direction
 * is a pure gradient, and therefore absorbable by the pressure without accelerating anything, only
 * if B has no horizontal variation. Where it does, the force has a curl, and a curl cannot be
 * balanced by any pressure field: it must spin up a circulation. The single exception is the one
 * every rotating planet uses — the Coriolis force acting on an existing zonal wind whose vertical
 * shear matches the horizontal density gradient. Cross-differentiating the meridional balance
 * against the hydrostatic one gives the condition:
 *
 *     dw/dr = (g / (f * rho_ref * r)) * drho/dtheta,        f = 2 Omega cos(theta)
 *
 * and its companion for the meridional wind, driven by the zonal density gradient:
 *
 *     dv/dr = -(g / (f * rho_ref * r sin(theta))) * drho/dphi
 *
 * ATJUP satisfies neither. Measured on the reference run, the model's own zonal wind has
 * essentially NO vertical shear — 0.0 to 0.26 m/s across the whole shell at every latitude — while
 * its temperature field demands 3 to 25 m/s. The jets are purely barotropic and the temperature
 * field is baroclinic, and the two were never required to know about each other, because until the
 * body forces were switched on nothing made them.
 *
 * WHAT IS ADJUSTED, AND WHY THAT ONE. The wind, not the temperature. Both would satisfy the
 * relation, but the temperature field carries the cloud decks and the whole condensation structure
 * that SaturationAdjustmentJup and PrecipitationJup are built on, and rebuilding it from the jets
 * gives a pole-to-equator contrast of order 150 K, which is not Jupiter. Adjusting the wind costs
 * 3 to 25 m/s of shear against jets of 130 m/s: a correction, not a replacement.
 *
 * The BAROTROPIC part of the jet is preserved exactly. Only the shear is replaced: the
 * mass-weighted column mean of the model's own zonal wind is kept, and the thermal-wind profile is
 * added as a deviation about its own column mean. So VelocityInitializerJup still sets how strong
 * each jet is and where it sits; this sets only how it varies with height.
 *
 * THE EQUATOR IS EXCLUDED, and not as a numerical convenience. f = 2 Omega cos(theta) vanishes
 * there, so the relation is singular — but more to the point geostrophy genuinely does not hold at
 * the equator, which is why Jupiter's equatorial jet is not a thermal-wind feature. Inside
 * ATJUP_TW_LAT_TAPER degrees of latitude the adjustment is blended smoothly back to the model's
 * own field, and the blend weight is reported.
 */

#pragma once

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

class cJupiterModel;

class ThermalWindJup {
public:
    explicit ThermalWindJup(cJupiterModel& model) : m(model) {}

    void run();

private:
    cJupiterModel& m;
};
