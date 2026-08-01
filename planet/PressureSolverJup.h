/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * Dynamic-pressure Poisson solve — ATJUP's binding of the SHARED implementation.
 * What used to differ between the two copies now lives either in PressureSolver.h or in
 * cJupiterModel: has_obstacle(), press_rigid_lid(), metricRadius() and
 * prepareProjectionBoundaries().
 */

#pragma once

#include "PressureSolver.h"

class cJupiterModel;

typedef PressureSolver<cJupiterModel> PressureSolverJup;
