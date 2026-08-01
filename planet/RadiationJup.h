/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * Grey multi-layer radiation — ATJUP's binding of the SHARED implementation in Radiation.h.
 * Jupiter's own radiative constants live in cJupiterModel.h, next to the rest of its parameters.
 */

#pragma once

#include "Radiation.h"

class cJupiterModel;

typedef Radiation<cJupiterModel> RadiationJup;
