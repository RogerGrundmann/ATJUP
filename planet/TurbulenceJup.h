/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * Turbulence closure — ATJUP's binding of the SHARED implementation in Turbulence.h.
 */

#pragma once

#include "Turbulence.h"

class cJupiterModel;

typedef Turbulence<cJupiterModel> TurbulenceJup;
