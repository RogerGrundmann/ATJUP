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
    // -----------------------------------------------------------------------
    static double clausius_clapeyron(double T_K, double A, double B){
        return std::exp(A / T_K + B);
    }

    static double saturation_vapour_pressure(double T_K,
            double C, double L0, double R, double del_alf, double del_bet){
        return std::exp(C
            + (-L0 / T_K + del_alf * std::log(T_K) + del_bet * T_K)
            / (1e-3 * R));
    }

    static double humility_critical(double x, double Hu_cr_max, double Hu_cr_mid){
        return (Hu_cr_max - Hu_cr_mid) * (x * x - 2.0 * x) + Hu_cr_max;
    }


private:
    cJupiterModel& m;

    static constexpr int iter_prec_end = 15;
//    static constexpr double q_diff_min = 1.0e-4;
    static constexpr double q_diff_min = 1.0e-3;
//    static constexpr double q_diff_min = 1.0e-2;
};

// Implementation is in SaturationAdjustmentJup.cpp
