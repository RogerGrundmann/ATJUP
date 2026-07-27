/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * Multi-layer grey-body radiation — STEP 1 SCAFFOLD.
 *
 * Declared as a friend of cJupiterModel so it may reach all private members
 * through the stored reference (same idiom as SaturationAdjustmentJup).
 *
 * Adapted from the ATOM_Precipitation MultiLayerRadiation.h, but the Earth
 * physics there (ice/snow surface-albedo feedback, ocean/land split, solar
 * surface-energy balance, Bignami water-vapour and Atwater-Ball CO2 laws,
 * and the surface-BC-entangled tridiagonal Thomas solve) does NOT apply to
 * Jupiter, so only the *architecture* is carried over:
 *   - per-column (j,k) grey layers, Stefan-Boltzmann emission sigma*T^4,
 *   - layer emissivity as an optical depth  eps_i = 1 - exp(-tau_i),
 *   - a two-stream up/down net-flux solve,
 *   - net-flux divergence -> radiative heating rate.
 *
 * Vertical orientation (verified in InitValues_Jup.cpp): i = i_base is the
 * deep bottom (~350 K, the intrinsic-heat boundary) and i = im-1 is the top
 * (radiates to space). i_base = i_topography[j][k] skips any solid GRS body.
 *
 * SCAFFOLD SCOPE (this step):
 *   - opacity is STUBBED to a constant per-layer optical depth (tau_stub); the
 *     real absorbers (H2-H2/H2-He CIA, CH4, NH3, clouds) come in later steps;
 *   - the result is written to DIAGNOSTIC arrays only (radiation, epsilon,
 *     Q_rad) and does NOT touch t or rhs_t — wiring the heating into the
 *     temperature equation is a later step, behind the same ATJUP_RADIATION knob.
 *   - In radiative equilibrium the interface net flux is uniform (= F_int) and
 *     Q_rad -> 0; a non-equilibrium column shows the tendency toward it. That is
 *     the intended self-test for this step.
*/

#pragma once

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

// cJupiterModel.h is included by the translation unit that uses this header.
class cJupiterModel;

class RadiationJup {
public:

    explicit RadiationJup(cJupiterModel& model) : m(model) {}

    // Physical constants / scaffold knobs.
    static constexpr double sigma  = 5.670374419e-8;  // Stefan-Boltzmann [W/m2/K4]
    static constexpr double F_int  = 5.4;             // Jupiter intrinsic heat flux [W/m2]

    // --- Absorbed shortwave (solar) ---
    // Jupiter emits ~13.9 W/m2 = 5.4 internal + ~8.5 absorbed solar. Injecting only F_int left
    // the column radiating ~2.6x what entered it, so it could never be in radiative balance —
    // matching the observed OLR was right for the wrong reason.
    //
    // Solar is handled as a SEPARATE SHORTWAVE CHANNEL, not by feeding the thermal Fd at the top.
    // The thermal opacity here is CIA-dominated (tau ~ P^2/T), which is the wrong absorber for
    // sunlight — sunlight is absorbed by CH4 bands and haze high up — and mixing the two would
    // also corrupt the net-flux diagnostic (Fu - Fd would net solar against thermal emission).
    //
    // Diurnally averaged insolation for a fast rotator with Jupiter's near-zero obliquity (3.13
    // deg) is S*(1-A)*cos(lat)/pi, whose area-weighted mean is exactly S*(1-A)/4 = 8.26 W/m2:
    //   mean = S(1-A)/pi * integral(cos^2) / integral(cos) = S(1-A)/pi * (pi/2)/2 = S(1-A)/4.
    // Equator gets 10.5 W/m2, poles 0.
    static constexpr double S_jup       = 50.5;    // solar constant at 5.2 AU [W/m2]
    static constexpr double albedo_bond = 0.343;   // Jupiter Bond albedo (so S*(1-A) is ABSORBED)
    // Shortwave optical depth per bar, measured down from the top. tau_sw = 1 near 1 bar puts
    // most of the absorption in the upper troposphere/stratosphere where CH4 and haze absorb.
    static constexpr double k_sw        = 1.0;     // [1/bar]

    // --- H2/He collision-induced absorption (CIA), grey parametrization (step 2) ---
    // CIA is the dominant thermal-IR opacity on the giant planets. It needs two collision
    // partners, each with number density n ∝ P/T, so the volume absorption ∝ (P/T)^2; over a
    // hydrostatic mass path dm = dP/g this gives a LAYER optical depth
    //     dtau = C_cia * comp * (P/T) * dP / g
    // i.e. the characteristic tau ∝ P^2/T of CIA. C_cia is a CALIBRATED grey coefficient
    // (lumps the frequency-integrated Borysow binary coefficient / k_B / mean molecular mass);
    // frequency-resolved tables are a later refinement. It is tuned so the thermal-IR
    // photosphere (tau=1 from the top) sits near ~0.5 bar, matching Jupiter.
    // comp = x_H2^2 + he_ratio*x_H2*x_He weights the H2-H2 and H2-He collision-pair
    // probabilities. A runtime multiplier ATJUP_CIA_STRENGTH (default 1) scales C_cia.
    static constexpr double P0_phys  = 1.0e5;     // physical pressure [Pa] at the p_stat=1 (T=t_ref) level ≈ 1 bar
    static constexpr double x_H2     = 0.863;     // H2 mole fraction
    static constexpr double x_He     = 0.134;     // He mole fraction
    static constexpr double he_ratio = 0.6;       // H2-He / H2-H2 grey binary-coefficient ratio
    static constexpr double C_cia    = 3.5e-6;    // calibrated grey CIA coefficient (photosphere ~0.5 bar)

    // --- CH4/NH3 gas bands + cloud/ice opacity, additive to the CIA tau (step 3) ---
    // Each absorber adds a layer optical depth  dtau = kappa * q * dP / g, where q is the
    // model's mass mixing ratio [kg/kg] and dP/g is the hydrostatic mass path [kg/m2]:
    //   - CH4 (7.7 um band) and NH3 (rotational + v2) grey band mass-absorption coefficients;
    //   - suspended cloud liquid (H2O/NH3/CH4 cloud) and ice (H2O/NH3/CH4 ice) grey continuum,
    //     after Stephens-type mass-absorption values.
    // The model over-condenses (as ATOM did: column condensate >> observed), which would make
    // every cloudy layer a blackbody and crush the OLR, so each layer's CLOUD tau is capped at
    // tau_cloud_cap. A runtime multiplier ATJUP_OPACITY_STRENGTH (default 1) scales the whole
    // non-CIA (gas+cloud) opacity for OLR tuning; CIA keeps its own ATJUP_CIA_STRENGTH.
    // NOTE: the model's CH4 field is large (q_ch4 ~ 0.246 kg/kg, ~100x the real Jovian value),
    // so kappa_ch4 is correspondingly SMALL to keep the CH4 optical depth physical — what is
    // calibrated is the grey OPTICAL DEPTH (kappa*q), not kappa alone. These give a combined
    // (CIA+gas+cloud) thermal-IR photosphere near ~0.25-0.35 bar and OLR ~ Jupiter's 14 W/m2.
    static constexpr double kappa_ch4   = 0.003;  // CH4 grey band mass opacity [m2/kg] (small: q_ch4 large)
    static constexpr double kappa_nh3   = 1.5;    // NH3 grey band mass opacity [m2/kg]
    static constexpr double kappa_cloud = 25.0;   // liquid-cloud grey mass opacity [m2/kg]
    static constexpr double kappa_ice   = 12.0;   // ice-cloud grey mass opacity [m2/kg]
    static constexpr double tau_cloud_cap = 0.4;  // per-layer cap on the cloud optical depth
    // OLR calibration factor on the whole gas+cloud opacity: with it the net radiative flux
    // at the emission level (~0.1-0.3 bar) is ~15 W/m2, matching Jupiter's ~14. (The true
    // top-of-atmosphere OLR is presently throttled below this by an anomalously cold model
    // top — an artifact of radiation not yet feeding T; revisit after the rhs_t coupling.)
    static constexpr double opac_cal    = 0.25;

    // Grey two-stream multi-layer solve. Fills radiation / epsilon / Q_rad.
    void run();

private:
    cJupiterModel& m;
};

// ---------------------------------------------------------------------------
// Implementation (header-only, like the ATOM original). Defined out-of-line so
// it can dereference the full cJupiterModel definition supplied by the includer.
// ---------------------------------------------------------------------------
#include "cJupiterModel.h"

inline void RadiationJup::run(){
    std::cout << std::endl << "      ATJUP: RadiationJup (grey CIA + CH4/NH3 + clouds)" << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    const int im = m.im, jm = m.jm, km = m.km;
    const double t_ref = m.t_ref;
    const double g     = m.g;

    // Composition-weighted collision-pair factor and the runtime CIA strength multiplier.
    const double comp = x_H2 * x_H2 + he_ratio * x_H2 * x_He;
    static const double cia_mult  = [](){ const char* e = getenv("ATJUP_CIA_STRENGTH");     return e ? atof(e) : 1.0; }();
    static const double opac_mult = [](){ const char* e = getenv("ATJUP_OPACITY_STRENGTH"); return e ? atof(e) : 1.0; }();
    const double cia_coeff = C_cia * cia_mult * comp / g;

    // Shortwave knobs. ATJUP_SOLAR=0 restores the pure-internal (F_int only) scaffold for A/B;
    // ATJUP_SOLAR_STRENGTH scales the absorbed flux; ATJUP_SW_TAU_PER_BAR moves the absorption
    // level (larger = absorbed higher up).
    static const int    solar_on   = [](){ const char* e = getenv("ATJUP_SOLAR");           return e ? atoi(e) : 1;   }();
    static const double solar_mult = [](){ const char* e = getenv("ATJUP_SOLAR_STRENGTH");  return e ? atof(e) : 1.0; }();
    static const double sw_tau_bar = [](){ const char* e = getenv("ATJUP_SW_TAU_PER_BAR");  return e ? atof(e) : k_sw; }();
    const double pi_ = 3.14159265358979323846;

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){

            const int i_base = std::max(m.i_topography[j][k], 0);  // deep boundary / first fluid cell
            const int i_top  = im - 1;                             // top boundary (to space)
            if(i_base >= i_top) continue;                          // degenerate column (all solid)

            // Thread-local column scratch.
            std::vector<double> B(im, 0.0);        // layer grey-body emission sigma*T^4 [W/m2]
            std::vector<double> eps(im, 0.0);      // layer emissivity
            std::vector<double> Fd(im + 1, 0.0);   // downward flux at bottom-interface of layer i (index i); top = im
            std::vector<double> Fu(im + 1, 0.0);   // upward   flux at bottom-interface of layer i (index i); top = im

            // Layer emission and H2/He CIA emissivity.
            for(int i = i_base; i <= i_top; i++){
                const double T = m.t.x[i][j][k] * t_ref;
                B[i] = sigma * T * T * T * T;

                // Physical pressures [Pa] from the model's dimensionless p_stat (~bars).
                const double P_lo = m.p_stat.x[i][j][k] * P0_phys;                     // layer bottom (higher P)
                const double P_hi = (i < i_top) ? m.p_stat.x[i + 1][j][k] * P0_phys    // layer top (lower P)
                                                : 0.0;                                 // top layer: vacuum above
                double dP = P_lo - P_hi;                                               // layer pressure thickness
                if(dP < 0.0) dP = 0.0;
                const double P_mean = 0.5 * (P_lo + P_hi);
                const double dm = dP / g;                                               // hydrostatic mass path [kg/m2]

                // CIA layer optical depth  tau_cia = cia_coeff * (P/T) * dP   (tau ∝ P^2/T).
                const double tau_cia = (T > 0.0) ? cia_coeff * (P_mean / T) * dP : 0.0;

                // CH4/NH3 gas band optical depth = kappa[m2/kg] * q * dm[kg/m2], which needs q
                // as a DIMENSIONLESS mass mixing ratio. The species fields are mass DENSITIES in
                // kg/m3 (see the units note at the top of printMinMax), so they are divided by
                // the local mixture density here. Before this they were fed in raw, i.e. a
                // density where a mixing ratio belongs, which understated tau wherever rho < 1 —
                // by a factor ~3 at the cloud decks and ~170 near the top of the shell.
                //
                // NOTE FOR RECALIBRATION: opac_cal = 0.25 was tuned against the old, wrong
                // quantity. ATJUP_OPACITY_STRENGTH is the lever if the photosphere level needs
                // to be brought back.
                const double rho_c = m.rho_mix.x[i][j][k];
                const double inv_rho = (rho_c > 0.0 && std::isfinite(rho_c)) ? 1.0 / rho_c : 0.0;
                const double q_ch4 = std::max(0.0, m.ch4.x[i][j][k]) * inv_rho;
                const double q_nh3 = std::max(0.0, m.nh3.x[i][j][k]) * inv_rho;
                const double tau_gas = opac_mult * opac_cal * (kappa_ch4 * q_ch4 + kappa_nh3 * q_nh3) * dm;

                // Cloud/ice continuum optical depth (all three condensing species), capped.
                // Same conversion as the gas bands: densities -> mixing ratios.
                const double q_liq = (std::max(0.0, m.h2o_cloud.x[i][j][k])
                                    + std::max(0.0, m.nh3_cloud.x[i][j][k])
                                    + std::max(0.0, m.ch4_cloud.x[i][j][k])) * inv_rho;
                const double q_ice = (std::max(0.0, m.h2o_ice.x[i][j][k])
                                    + std::max(0.0, m.nh3_ice.x[i][j][k])
                                    + std::max(0.0, m.ch4_ice.x[i][j][k])) * inv_rho;
                double tau_cloud = opac_mult * opac_cal * (kappa_cloud * q_liq + kappa_ice * q_ice) * dm;
                if(tau_cloud > tau_cloud_cap) tau_cloud = tau_cloud_cap;

                const double tau = tau_cia + tau_gas + tau_cloud;
                eps[i] = 1.0 - std::exp(-tau);
            }

            // --- Shortwave (solar) channel, independent of the thermal sweeps. ---
            // Absorbed insolation at the top of this column, then Beer-Lambert attenuation with
            // an optical depth measured DOWN FROM THE TOP so that Fsw at the top interface is
            // exactly F_sun_toa. Energy is then conserved by construction: everything absorbed in
            // the layers plus the residual reaching the deep boundary sums to F_sun_toa exactly.
            std::vector<double> Fsw(im + 2, 0.0);   // downward shortwave at bottom-interface of layer i
            double F_sun_toa = 0.0;
            if(solar_on != 0){
                // cos(latitude) = sin(colatitude); colatitude = j*pi/(jm-1), so j=0/jm-1 are the
                // poles (cos_lat = 0) and j=(jm-1)/2 the equator (cos_lat = 1).
                const double colat   = pi_ * (double)j / (double)(jm - 1);
                const double cos_lat = std::max(0.0, std::sin(colat));
                F_sun_toa = solar_mult * S_jup * (1.0 - albedo_bond) * cos_lat / pi_;

                // p_stat is ~bars; the interface above layer i_top is vacuum (p = 0), which is
                // also the zero point of the shortwave optical depth.
                for(int i = i_base; i <= i_top + 1; i++){
                    const double p_bar = (i <= i_top) ? m.p_stat.x[i][j][k] : 0.0;
                    Fsw[i] = F_sun_toa * std::exp(-sw_tau_bar * std::max(0.0, p_bar));
                }
            }

            // Downward sweep (top -> bottom). Fd[i] is the downward flux leaving the
            // bottom of layer i; Fd[i_top+1] is the incoming THERMAL flux at the top
            // boundary, which is zero — space is cold, and the solar term is the separate
            // shortwave channel above, not a downward thermal flux.
            Fd[i_top + 1] = 0.0;
            for(int i = i_top; i >= i_base; i--)
                Fd[i] = Fd[i + 1] * (1.0 - eps[i]) + eps[i] * B[i];

            // Bottom boundary: net flux up through the deep interface must carry the intrinsic
            // flux PLUS whatever shortwave survived to depth (absorbed there and re-emitted in
            // the thermal channel), so the deep boundary emits Fd_bottom + F_int + Fsw_bottom.
            Fu[i_base] = Fd[i_base] + F_int + Fsw[i_base];

            // Upward sweep (bottom -> top). Fu[i+1] is the upward flux leaving the
            // top of layer i (= entering the bottom of layer i+1).
            for(int i = i_base; i <= i_top; i++)
                Fu[i + 1] = Fu[i] * (1.0 - eps[i]) + eps[i] * B[i];

            // Interface net flux (up - down) and per-layer heating from its divergence.
            // In radiative equilibrium every interface net flux equals F_int and the
            // divergence (hence Q_rad) vanishes.
            for(int i = i_base; i <= i_top; i++){
                const double net_bot = Fu[i]     - Fd[i];       // net flux at bottom interface of layer i
                const double net_top = Fu[i + 1] - Fd[i + 1];   // net flux at top interface of layer i

                // Layer thickness in METRES. get_layer_height() is in km (L_atm is specified in
                // km), so the raw difference would make Q_rad 1000x too large.
                double dz = (i < i_top) ? m.layer_thickness_m(i) : m.layer_thickness_m(i - 1);
                if(dz <= 0.0) dz = 1.0;   // guard degenerate layer heights

                // Shortwave heating: convergence of the downward solar beam in this layer.
                // Fsw decreases downward, so (Fsw[i+1] - Fsw[i]) >= 0 is absorbed energy and the
                // contribution to Q_rad is a HEATING term.
                const double sw_absorbed = Fsw[i + 1] - Fsw[i];          // [W/m2] absorbed in layer i
                const double Q_sw        = sw_absorbed / dz;             // [W/m3]

                m.epsilon.x[i][j][k]   = eps[i];
                m.radiation.x[i][j][k] = 0.5 * (net_bot + net_top);      // layer-centre net THERMAL flux [W/m2]
                m.Q_rad.x[i][j][k]     = (net_bot - net_top) / dz + Q_sw;  // total radiative heating [W/m3]
            }

            // Copy the deepest fluid value into any solid cells below so plots/BCs
            // see a filled column.
            for(int i = 0; i < i_base; i++){
                m.epsilon.x[i][j][k]   = m.epsilon.x[i_base][j][k];
                m.radiation.x[i][j][k] = m.radiation.x[i_base][j][k];
                m.Q_rad.x[i][j][k]     = 0.0;
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for RadiationJup\n", elapsed.count() * 1e-9);
    std::cout << "      ATJUP: RadiationJup ended" << std::endl;
}
