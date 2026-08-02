/*
 * Atmosphere General Circulation Modell (ATJUP)
 * Standalone saturation-adjustment and microphysics class for the Jupiter model.
 * Declared as friend of cJupiterModel so it may access all private members
 * through the stored reference.
 *
 * Reference algorithm:
 *   Tao, W.-K., Simpson, J., and McCumber, M.:
 *   "An Ice-Water Saturation Adjustment", AMS Notes and Correspondence, 1988.
*/

#pragma once

#include "ATPhys.h"   // saturation_vapour_pressure, clausius_clapeyron

#include <cmath>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <string>
#include <iostream>

// cJupiterModel.h is included by the translation unit that uses this header.
// We only need the forward declaration here to avoid a circular include.
class cJupiterModel;
// Array is used by value (pointer indirection only), so forward-declare it too.
class Array;

using namespace std;

class SaturationAdjustmentJup {
public:

    explicit SaturationAdjustmentJup(cJupiterModel& model) : m(model) {}


    // -----------------------------------------------------------------------
    // Saturation adjustment (Tao et al. 1988)
    // Adjusts vapour / cloud / ice fields to thermodynamic equilibrium and
    // updates t and p_stat in-place.
    // -----------------------------------------------------------------------
    void run(const std::string& gas,
             double coeff_A,   double coeff_B,
             double coeff_A_i, double coeff_B_i,
             double t_0,       double t_00,
             double ep,        double lv,  double ls,
             double cp_gas,    double r_gas,
             double C,         double L0,  double R,
             double del_alf,   double del_bet,   double m_mol,
             double C_i,       double L0_i,
             double del_alf_i, double del_bet_i,
             Array& c,         Array& cloud,   Array& ice);


    // -----------------------------------------------------------------------
    // Static helper functions — usable without a SaturationAdjustmentJup instance
    // The three static helpers that used to sit here are gone. saturation_vapour_pressure and
    // clausius_clapeyron are ATPhys:: functions now — one implementation for every planet instead
    // of a copy per class — and humility_critical had no caller in either model, which is why it
    // was deleted rather than moved: a shared header is the wrong place to preserve dead code.
    // (cJupiterModel::Humility_critical, capital H, is a DIFFERENT function.)


private:
    cJupiterModel& m;

    static constexpr int iter_prec_end = 15;
//    static constexpr double q_diff_min = 1.0e-4;
    static constexpr double q_diff_min = 1.0e-3;
//    static constexpr double q_diff_min = 1.0e-2;
};

// Implementation is in SaturationAdjustmentJup.cpp
