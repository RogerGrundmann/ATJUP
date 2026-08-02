/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * Saturation adjustment and mixed-phase partitioning — ATJUP's binding of the SHARED
 * implementation. The algorithm, and the long note on why it is shared, are in
 * SaturationAdjustment.h.
 *
 * What used to differ between the two copies now lives either in that file or in cJupiterModel:
 * is_solid() for the obstacle interior, and satadj_updates_pstat() for whether the hydrostatic
 * pressure is rebuilt from the adjusted temperature (ATJUP: no).
 *
 * Reference algorithm:
 *   Tao, W.-K., Simpson, J., and McCumber, M.:
 *   "An Ice-Water Saturation Adjustment", AMS Notes and Correspondence, 1988.
 */

#pragma once

#include "SaturationAdjustment.h"

class cJupiterModel;

typedef SaturationAdjustment<cJupiterModel> SaturationAdjustmentJup;
