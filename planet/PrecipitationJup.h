/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * Ice / precipitation microphysics — ATJUP's binding of the SHARED implementation.
 */

#pragma once

#include "Precipitation.h"

class cJupiterModel;

typedef Precipitation<cJupiterModel> PrecipitationJup;
