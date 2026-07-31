/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * Dry convective adjustment — ATJUP's binding of the SHARED implementation.
 *
 * The algorithm, the reasoning and the knobs live in ConvectiveAdjustment.h, which is
 * byte-identical in every model that uses it and knows nothing about which planet it runs on.
 * This file exists only so the call sites keep their familiar name.
 */

#pragma once

#include "ConvectiveAdjustment.h"

class cJupiterModel;

typedef ConvectiveAdjustment<cJupiterModel> ConvectiveAdjustmentJup;
