/*
 * Atmosphere General Circulation Modell(ATJUP) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include "cJupiterModel.h"

#include <cstdlib>   // getenv/atof for the radiative-coupling knob

using namespace std;


void cJupiterModel::RHSJup(int i, int j, int k, const CellGeometry& geo){

    // All geometric quantities come from the precomputed struct —
    // NO sin(), cos(), division, or reciprocal computation here.
    const double rm                  = geo.rm;
    const double sinthe              = geo.sinthe;
    const double costhe              = geo.costhe;
    const double cotanthe            = geo.cotanthe;
    const double inv_rm              = geo.inv_rm;
    const double inv_rm2             = geo.inv_rm2;
    const double inv_rmsinthe        = geo.inv_rmsinthe;
    const double inv_rm2sinthe       = geo.inv_rm2sinthe;
    const double inv_rm2sinthe2      = geo.inv_rm2sinthe2;
    const double costhe_inv_rm2sinthe = geo.costhe_inv_rm2sinthe;

    const double inv_2dr   = geo.inv_2dr;
    const double inv_2dthe = geo.inv_2dthe;
    const double inv_2dphi = geo.inv_2dphi;
    const double inv_dr2   = geo.inv_dr2;
    const double inv_dthe2 = geo.inv_dthe2;
    const double inv_dphi2 = geo.inv_dphi2;
    const double exp_rm    = geo.exp_rm;
    const double exp_2_rm  = geo.exp_2_rm;

    // Cache local cell values
    const double u_ijk = u.x[i][j][k];
    const double v_ijk = v.x[i][j][k];
    const double w_ijk = w.x[i][j][k];

    // ---- First-order derivative storage ----
    double dudr, dvdr, dwdr, dtdr, dpdr;
    double dh2odr, dh2ocdr, dh2oidr;
    double dh2sdr;
    double dnh3dr, dnh3cdr, dnh3idr, dnh4shdr;
    double dch4dr, dch4cdr, dch4idr;

    double dudthe, dvdthe, dwdthe, dtdthe, dpdthe;
    double dh2odthe, dh2ocdthe, dh2oidthe;
    double dh2sdthe;
    double dnh3dthe, dnh3cdthe, dnh3idthe, dnh4shdthe;
    double dch4dthe, dch4cdthe, dch4idthe;

    double dudphi, dvdphi, dwdphi, dtdphi, dpdphi;
    double dh2odphi, dh2ocdphi, dh2oidphi;
    double dh2sdphi;
    double dnh3dphi, dnh3cdphi, dnh3idphi, dnh4shdphi;
    double dch4dphi, dch4cdphi, dch4idphi;

    // ---- Second-order derivative storage ----
    double d2udr2, d2vdr2, d2wdr2, d2tdr2;
    double d2h2odr2, d2h2ocdr2, d2h2oidr2;
    double d2h2sdr2;
    double d2nh3dr2, d2nh3cdr2, d2nh3idr2, d2nh4shdr2;
    double d2ch4dr2, d2ch4cdr2, d2ch4idr2;

    double d2udthe2, d2vdthe2, d2wdthe2, d2tdthe2;
    double d2h2odthe2, d2h2ocdthe2, d2h2oidthe2;
    double d2h2sdthe2;
    double d2nh3dthe2, d2nh3cdthe2, d2nh3idthe2, d2nh4shdthe2;
    double d2ch4dthe2, d2ch4cdthe2, d2ch4idthe2;

    double d2udphi2, d2vdphi2, d2wdphi2, d2tdphi2;
    double d2h2odphi2, d2h2ocdphi2, d2h2oidphi2;
    double d2h2sdphi2;
    double d2nh3dphi2, d2nh3cdphi2, d2nh3idphi2, d2nh4shdphi2;
    double d2ch4dphi2, d2ch4cdphi2, d2ch4idphi2;


    // ===== R-direction derivatives (central differences) =====
    // exp_rm  = 1/(rm+1)   scales the first  derivative: ∂f/∂r_phys ≈ (∂f/∂r_model)*exp_rm
    // exp_2_rm = exp_rm²   scales the second derivative: ∂²f/∂r_phys² ≈ (∂²f/∂r_model²)*exp_2_rm
    #define COMPUTE_DR(FIELD, d1, d2) \
        d1 = (FIELD.x[i+1][j][k] - FIELD.x[i-1][j][k]) * inv_2dr * exp_rm; \
        d2 = (FIELD.x[i+1][j][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i-1][j][k]) * inv_dr2 * exp_2_rm;

    COMPUTE_DR(u,         dudr,    d2udr2)
    COMPUTE_DR(v,         dvdr,    d2vdr2)
    COMPUTE_DR(w,         dwdr,    d2wdr2)
    COMPUTE_DR(t,         dtdr,    d2tdr2)
    COMPUTE_DR(h2o,       dh2odr,  d2h2odr2)
    COMPUTE_DR(h2o_cloud, dh2ocdr, d2h2ocdr2)
    COMPUTE_DR(h2o_ice,   dh2oidr, d2h2oidr2)
    COMPUTE_DR(h2s,       dh2sdr,  d2h2sdr2)
    COMPUTE_DR(nh3,       dnh3dr,  d2nh3dr2)
    COMPUTE_DR(nh3_cloud, dnh3cdr, d2nh3cdr2)
    COMPUTE_DR(nh3_ice,   dnh3idr, d2nh3idr2)
    COMPUTE_DR(ch4,       dch4dr,  d2ch4dr2)
    COMPUTE_DR(ch4_cloud, dch4cdr, d2ch4cdr2)
    COMPUTE_DR(ch4_ice,   dch4idr, d2ch4idr2)
    COMPUTE_DR(nh4sh,     dnh4shdr,d2nh4shdr2)
    dpdr = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) * inv_2dr * exp_rm;
    #undef COMPUTE_DR


    // ===== Theta-direction derivatives (central differences) =====
    #define COMPUTE_DTHE(FIELD, d1, d2) \
        d1 = (FIELD.x[i][j+1][k] - FIELD.x[i][j-1][k]) * inv_2dthe; \
        d2 = (FIELD.x[i][j+1][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j-1][k]) * inv_dthe2;

    COMPUTE_DTHE(u,         dudthe,    d2udthe2)
    COMPUTE_DTHE(v,         dvdthe,    d2vdthe2)
    COMPUTE_DTHE(w,         dwdthe,    d2wdthe2)
    COMPUTE_DTHE(t,         dtdthe,    d2tdthe2)
    COMPUTE_DTHE(h2o,       dh2odthe,  d2h2odthe2)
    COMPUTE_DTHE(h2o_cloud, dh2ocdthe, d2h2ocdthe2)
    COMPUTE_DTHE(h2o_ice,   dh2oidthe, d2h2oidthe2)
    COMPUTE_DTHE(h2s,       dh2sdthe,  d2h2sdthe2)
    COMPUTE_DTHE(nh3,       dnh3dthe,  d2nh3dthe2)
    COMPUTE_DTHE(nh3_cloud, dnh3cdthe, d2nh3cdthe2)
    COMPUTE_DTHE(nh3_ice,   dnh3idthe, d2nh3idthe2)
    COMPUTE_DTHE(ch4,       dch4dthe,  d2ch4dthe2)
    COMPUTE_DTHE(ch4_cloud, dch4cdthe, d2ch4cdthe2)
    COMPUTE_DTHE(ch4_ice,   dch4idthe, d2ch4idthe2)
    COMPUTE_DTHE(nh4sh,     dnh4shdthe,d2nh4shdthe2)
    dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) * inv_2dthe;
    #undef COMPUTE_DTHE


    // ===== Phi-direction derivatives (central differences) =====
    #define COMPUTE_DPHI(FIELD, d1, d2) \
        d1 = (FIELD.x[i][j][k+1] - FIELD.x[i][j][k-1]) * inv_2dphi; \
        d2 = (FIELD.x[i][j][k+1] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j][k-1]) * inv_dphi2;

    COMPUTE_DPHI(u,         dudphi,    d2udphi2)
    COMPUTE_DPHI(v,         dvdphi,    d2vdphi2)
    COMPUTE_DPHI(w,         dwdphi,    d2wdphi2)
    COMPUTE_DPHI(t,         dtdphi,    d2tdphi2)
    COMPUTE_DPHI(h2o,       dh2odphi,  d2h2odphi2)
    COMPUTE_DPHI(h2o_cloud, dh2ocdphi, d2h2ocdphi2)
    COMPUTE_DPHI(h2o_ice,   dh2oidphi, d2h2oidphi2)
    COMPUTE_DPHI(h2s,       dh2sdphi,  d2h2sdphi2)
    COMPUTE_DPHI(nh3,       dnh3dphi,  d2nh3dphi2)
    COMPUTE_DPHI(nh3_cloud, dnh3cdphi, d2nh3cdphi2)
    COMPUTE_DPHI(nh3_ice,   dnh3idphi, d2nh3idphi2)
    COMPUTE_DPHI(ch4,       dch4dphi,  d2ch4dphi2)
    COMPUTE_DPHI(ch4_cloud, dch4cdphi, d2ch4cdphi2)
    COMPUTE_DPHI(ch4_ice,   dch4idphi, d2ch4idphi2)
    COMPUTE_DPHI(nh4sh,     dnh4shdphi,d2nh4shdphi2)
    dpdphi = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dphi;
    #undef COMPUTE_DPHI


    // ===== Coriolis and centrifugal forces =====
    double Coriolis_rad  = -2.0 * omega * sinthe * w_ijk;
    double Coriolis_the  = +2.0 * omega * costhe * w_ijk;
    double Coriolis_phi  = +2.0 * omega * (-costhe * v_ijk + sinthe * u_ijk);

    double centrifugal_rad = omega * omega * rm;
    double centrifugal_the = omega * omega * rm * fabs(sinthe);

    double coeff_energy_p = u_0 * u_0 / (cp_mix * t_ref);


    // ===== Transport terms (advection) =====
    double v_invrm = v_ijk * inv_rm;
    double w_invrs = w_ijk * inv_rmsinthe;

    double pressure_t = coeff_energy_p
        * (u_ijk * dpdr + v_invrm * dpdthe + w_invrs * dpdphi);

    double transport_t = u_ijk * dtdr + v_invrm * dtdthe + w_invrs * dtdphi;

    double transport_u = u_ijk * dudr + v_invrm * dudthe + w_invrs * dudphi
        - (v_ijk * v_ijk + w_ijk * w_ijk) * inv_rm;
    double transport_v = u_ijk * dvdr + v_invrm * dvdthe + w_invrs * dvdphi
        + (u_ijk * v_ijk - w_ijk * w_ijk * cotanthe) * inv_rm;
    double transport_w = u_ijk * dwdr + v_invrm * dwdthe + w_invrs * dwdphi
        + (w_ijk * u_ijk + v_ijk * w_ijk * cotanthe) * inv_rm;

    double transport_h2o       = u_ijk * dh2odr  + v_invrm * dh2odthe  + w_invrs * dh2odphi;
    double transport_h2o_cloud = u_ijk * dh2ocdr + v_invrm * dh2ocdthe + w_invrs * dh2ocdphi;
    double transport_h2o_ice   = u_ijk * dh2oidr + v_invrm * dh2oidthe + w_invrs * dh2oidphi;

    double transport_h2s       = u_ijk * dh2sdr  + v_invrm * dh2sdthe  + w_invrs * dh2sdphi;

    double transport_nh3       = u_ijk * dnh3dr  + v_invrm * dnh3dthe  + w_invrs * dnh3dphi;
    double transport_nh3_cloud = u_ijk * dnh3cdr + v_invrm * dnh3cdthe + w_invrs * dnh3cdphi;
    double transport_nh3_ice   = u_ijk * dnh3idr + v_invrm * dnh3idthe + w_invrs * dnh3idphi;

    double transport_ch4       = u_ijk * dch4dr  + v_invrm * dch4dthe  + w_invrs * dch4dphi;
    double transport_ch4_cloud = u_ijk * dch4cdr + v_invrm * dch4cdthe + w_invrs * dch4cdphi;
    double transport_ch4_ice   = u_ijk * dch4idr + v_invrm * dch4idthe + w_invrs * dch4idphi;

    double transport_nh4sh     = u_ijk * dnh4shdr + v_invrm * dnh4shdthe + w_invrs * dnh4shdphi;


    // ===== Diffusion terms =====
    double two_inv_rm    = 2.0 * inv_rm;
    double v_metric      = (1.0 + costhe / (geo.sinthe2)) * inv_rm2;

    double diffusion_t = d2tdr2       + dtdr       * two_inv_rm + d2tdthe2       * inv_rm2
        + dtdthe       * costhe_inv_rm2sinthe + d2tdphi2       * inv_rm2sinthe2;

    double diffusion_u = d2udr2 + 2.0 * u_ijk * inv_rm2 + d2udthe2 * inv_rm2
        + 4.0 * dudr * inv_rm  + dudthe * costhe_inv_rm2sinthe
        + d2udphi2 * inv_rm2sinthe2;
    double diffusion_v = d2vdr2       + dvdr       * two_inv_rm + d2vdthe2       * inv_rm2
        + dvdthe       * costhe_inv_rm2sinthe - v_metric * v_ijk
        + d2vdphi2     * inv_rm2sinthe2
        + 2.0 * dudthe * inv_rm2
        - dwdphi * 2.0 * costhe * inv_rm2sinthe2;
    double diffusion_w = d2wdr2       + dwdr       * two_inv_rm + d2wdthe2       * inv_rm2
        + dwdthe       * costhe_inv_rm2sinthe - v_metric * w_ijk
        + d2wdphi2     * inv_rm2sinthe2
        + 2.0 * dudphi * inv_rm2sinthe
        + dvdphi * 2.0 * costhe * inv_rm2sinthe2;

    double diffusion_h2o = d2h2odr2  + dh2odr  * two_inv_rm + d2h2odthe2  * inv_rm2
        + dh2odthe  * costhe_inv_rm2sinthe + d2h2odphi2  * inv_rm2sinthe2;
    double diffusion_h2o_cloud = d2h2ocdr2 + dh2ocdr * two_inv_rm + d2h2ocdthe2 * inv_rm2
        + dh2ocdthe * costhe_inv_rm2sinthe + d2h2ocdphi2 * inv_rm2sinthe2;
    double diffusion_h2o_ice = d2h2oidr2 + dh2oidr * two_inv_rm + d2h2oidthe2 * inv_rm2
        + dh2oidthe * costhe_inv_rm2sinthe + d2h2oidphi2 * inv_rm2sinthe2;

    double diffusion_h2s = d2h2sdr2  + dh2sdr  * two_inv_rm + d2h2sdthe2  * inv_rm2
        + dh2sdthe  * costhe_inv_rm2sinthe + d2h2sdphi2  * inv_rm2sinthe2;

    double diffusion_nh3 = d2nh3dr2  + dnh3dr  * two_inv_rm + d2nh3dthe2  * inv_rm2
        + dnh3dthe  * costhe_inv_rm2sinthe + d2nh3dphi2  * inv_rm2sinthe2;
    double diffusion_nh3_cloud = d2nh3cdr2 + dnh3cdr * two_inv_rm + d2nh3cdthe2 * inv_rm2
        + dnh3cdthe * costhe_inv_rm2sinthe + d2nh3cdphi2 * inv_rm2sinthe2;
    double diffusion_nh3_ice = d2nh3idr2 + dnh3idr * two_inv_rm + d2nh3idthe2 * inv_rm2
        + dnh3idthe * costhe_inv_rm2sinthe + d2nh3idphi2 * inv_rm2sinthe2;

    double diffusion_ch4 = d2ch4dr2  + dch4dr  * two_inv_rm + d2ch4dthe2  * inv_rm2
        + dch4dthe  * costhe_inv_rm2sinthe + d2ch4dphi2  * inv_rm2sinthe2;
    double diffusion_ch4_cloud = d2ch4cdr2 + dch4cdr * two_inv_rm + d2ch4cdthe2 * inv_rm2
        + dch4cdthe * costhe_inv_rm2sinthe + d2ch4cdphi2 * inv_rm2sinthe2;
    double diffusion_ch4_ice = d2ch4idr2 + dch4idr * two_inv_rm + d2ch4idthe2 * inv_rm2
        + dch4idthe * costhe_inv_rm2sinthe + d2ch4idphi2 * inv_rm2sinthe2;

    double diffusion_nh4sh = d2nh4shdr2 + dnh4shdr * two_inv_rm + d2nh4shdthe2 * inv_rm2
        + dnh4shdthe * costhe_inv_rm2sinthe + d2nh4shdphi2 * inv_rm2sinthe2;


    // ===== RHS assembly =====
    double dpdr_term   = dpdr;
    double dpdthe_term = dpdthe * inv_rm;
    double dpdphi_term = dpdphi * inv_rmsinthe;

    // ===== Radiative heating source (step 4, opt-in) =====
    // Convert the diagnostic radiative flux divergence Q_rad [W/m3] (RadiationJup) into a
    // nondimensional temperature tendency and add it to rhs_t. Physically dT/dt = Q_rad/(rho*cp);
    // nondimensionalised by the energy-equation scaling (radial length L_rad, velocity u_0,
    // temperature t_ref):
    //   radiation_t = rad_coupling * Q_rad * L_rad / (rho * cp_mix * u_0 * t_ref).
    // rho is the LOCAL density from the ideal-gas law (essential so the thin, cold upper
    // atmosphere — where Q_rad>0 heats — responds strongly and relaxes toward radiative
    // equilibrium). Gated by ATJUP_RAD_COUPLING (default 0 = off, bit-identical); Q_rad is
    // nonzero only when ATJUP_RADIATION is enabled.
    //
    // SCALING: rad_coupling = 1.0 is the PHYSICALLY CORRECT value — the expression above is the
    // exact nondimensional form of dT/dt = Q/(rho*cp) under this model's scaling (lengths by
    // L_rad, velocity u_0, temperature t_ref). Verified: Q_rad~5e-3 W/m3 at the cloud decks
    // gives dT/dt = Q/(rho*cp) ~ 4.3e-6 K/s, and radiation_t*dt*t_ref reproduces that per step.
    // Values >> 1 do NOT correct a scaling error; they are a deliberate ACCELERATION factor.
    // The reason one is tempting: dt = 0.001 nondimensional is only dt*L_rad/u_0 ~ 1.4 s of
    // Jupiter time, so a 100-iteration run spans ~140 s while the radiative relaxation time is
    // ~1e7 s. At coupling=1 radiative equilibration therefore needs ~1e7/1.4 ~ 7e6 iterations;
    // raising the knob buys that equilibration in fewer steps at the cost of a distorted ratio
    // between the radiative and advective timescales. (An earlier +8 K/100-iter result at
    // coupling=10 was mostly a units bug: Q_rad was computed with a layer thickness in km
    // instead of m and so was 1000x too large, making the effective multiplier ~1e4.)
    static const double rad_coupling = [](){ const char* e = getenv("ATJUP_RAD_COUPLING"); return e ? atof(e) : 0.0; }();
    double radiation_t = 0.0;
    if(rad_coupling != 0.0){
        constexpr double R_H2He = 3600.0;                       // specific gas constant of the H2/He mix [J/(kg*K)]
        const double T_phys = t.x[i][j][k] * t_ref;             // [K]
        const double P_phys = p_stat.x[i][j][k] * 1.0e5;        // p_stat ~ bars -> [Pa]
        const double rho    = (T_phys > 1.0) ? P_phys / (R_H2He * T_phys) : 0.0;  // [kg/m3]
        const double L_rad  = L_atm * 1.0e3;                    // atmosphere thickness [m]
        if(rho > 0.0 && cp_mix > 0.0){
            radiation_t = rad_coupling * Q_rad.x[i][j][k] * L_rad
                        / (rho * cp_mix * u_0 * t_ref);
            // Explicit-scheme stability limiter: in a very low-density cell (thin upper
            // atmosphere / deep polar corner) the 1/rho factor can make the tendency blow up
            // and destabilise the pole. Guard non-finite first (the cap's >/< tests are both
            // false for NaN and would let it through), then cap — normal top-of-atmosphere
            // values are ~0.05, so this only bites on the runaway, preserving the heating sign.
            constexpr double rad_t_max = 0.5;
            if(!std::isfinite(radiation_t)) radiation_t = 0.0;
            else if(radiation_t >  rad_t_max) radiation_t =  rad_t_max;
            else if(radiation_t < -rad_t_max) radiation_t = -rad_t_max;
        }
    }

    // ===== Latent-heat source from precipitation microphysics (phase 2c, opt-in) =====
    // Same conversion as the radiative term: Q_precip [W/m3] (PrecipitationJup, summed over
    // H2O/NH3/NH4SH) -> nondimensional temperature tendency via dT/dt = Q/(rho*cp), scaled by
    // L_rad/(u_0*t_ref). Positive where riming/freezing release fusion heat, negative where
    // melting or rain evaporation absorb it. Gated by ATJUP_PRECIP_COUPLING (default 0 =
    // bit-identical); Q_precip is nonzero only when ATJUP_PRECIP is enabled.
    // NOTE the density used here is r_mix, NOT the local ideal-gas density used by the
    // radiative term. Q_rad comes from real radiative fluxes, so dividing it by the local
    // density is correct. Q_precip instead derives from the condensate fields, which
    // SaturationAdjustmentJup defines as mixing ratios scaled by the REFERENCE density r_mix
    // (and whose own latent heat it converts with /(cp_mix*r_mix)). The r_mix therefore
    // cancels, leaving the true mixing-ratio tendency; using the local density here would
    // instead inflate the heating by r_mix/rho_local (~13x at the cloud decks).
    static const double precip_coupling = [](){ const char* e = getenv("ATJUP_PRECIP_COUPLING"); return e ? atof(e) : 0.0; }();
    double precip_t = 0.0;
    if(precip_coupling != 0.0){
        const double L_rad = L_atm * 1.0e3;                     // atmosphere thickness [m]
        if(r_mix > 0.0 && cp_mix > 0.0){
            precip_t = precip_coupling * Q_precip.x[i][j][k] * L_rad
                     / (r_mix * cp_mix * u_0 * t_ref);
            // Same explicit-scheme guard as the radiative term. r_mix is a constant here so
            // there is no 1/rho blow-up, but latent heating is a stiff, locally concentrated
            // source (it switches on hard at the freezing level), so keep the limiter.
            constexpr double precip_t_max = 0.5;
            if(!std::isfinite(precip_t)) precip_t = 0.0;
            else if(precip_t >  precip_t_max) precip_t =  precip_t_max;
            else if(precip_t < -precip_t_max) precip_t = -precip_t_max;
        }
    }

    rhs_t.x[i][j][k] =
        + pressure_t
        - transport_t
        + diffusion_t / (re * pr)
        + radiation_t
        + precip_t;

    rhs_u.x[i][j][k] =
        - dpdr_term
//        + buoyancy * g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])
//                      / (r_mix * R_mix * t.x[i][j][k] * t_ref)
        + buoyancy * g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])
                      / (r_mix * R_mix * t.x[i][j][k] * t_ref)
        - transport_u
        + diffusion_u / re
        - Coriolis    * Coriolis_rad
        - centrifugal * centrifugal_rad;


/*
    // 1. Lokale Einstein-Viskosität berechnen (dimensionslos)
    // mu_base entspricht 1.0, da re bereits im Nenner des Gesamterms steht
    double mu_eff_cell = 1.0 * (1.0 + 2.5 * rho_mix[i][j][k]/rg_nh4sh);

    // 2. Diffusionsterm mit variabler Viskosität (harmonisch gemittelt)
    // Hier wird mu_eff an den Flächen i+1/2 und i-1/2 berechnet
    double mu_east = (2.0 * mu_eff[i][j][k] * mu_eff[i+1][j][k]) / (mu_eff[i][j][k] + mu_eff[i+1][j][k] + 1e-20);
    double mu_west = (2.0 * mu_eff[i][j][k] * mu_eff[i-1][j][k]) / (mu_eff[i][j][k] + mu_eff[i-1][j][k] + 1e-20);

    // Der neue Diffusionsterm (ersetzt dein altes diffusion_u / re)
    double variable_diffusion = (mu_east * (u[i+1] - u[i]) - mu_west * (u[i] - u[i-1])) / (re * dx * dx);

    // 3. RHS Zusammensetzung
    rhs_u.x[i][j][k] = 
        - dpdr_term 
        + buoyancy_term
        - transport_u 
        + variable_diffusion  // <--- Das ist die Änderung
        - Coriolis_term 
        - centrifugal_term;

*/


    rhs_v.x[i][j][k] =
        - dpdthe_term
        - transport_v
        + diffusion_v / re
        - Coriolis    * Coriolis_the
        - centrifugal * centrifugal_the;

    rhs_w.x[i][j][k] =
        - dpdphi_term
        - transport_w
        + diffusion_w / re
        - Coriolis    * Coriolis_phi;

    rhs_h2o.x[i][j][k] =
        - transport_h2o
        + diffusion_h2o / (sc_h2o * re);

    rhs_h2o_cloud.x[i][j][k] =
        - transport_h2o_cloud
        + diffusion_h2o_cloud / (sc_h2o * re);

    rhs_h2o_ice.x[i][j][k] =
        - transport_h2o_ice
        + diffusion_h2o_ice / (sc_h2o * re);

    rhs_h2s.x[i][j][k] =
        - transport_h2s
        + diffusion_h2s / (sc_h2s * re)
        + chemical_reaction * massflux_h2s.x[i][j][k];

    rhs_nh3.x[i][j][k] =
        - transport_nh3
        + diffusion_nh3 / (sc_nh3 * re)
        + chemical_reaction * massflux_nh3.x[i][j][k];

    rhs_nh3_cloud.x[i][j][k] =
        - transport_nh3_cloud
        + diffusion_nh3_cloud / (sc_nh3 * re);

    rhs_nh3_ice.x[i][j][k] =
        - transport_nh3_ice
        + diffusion_nh3_ice / (sc_nh3 * re);

    rhs_ch4.x[i][j][k] =
        - transport_ch4
        + diffusion_ch4 / (sc_ch4 * re);

    rhs_ch4_cloud.x[i][j][k] =
        - transport_ch4_cloud
        + diffusion_ch4_cloud / (sc_ch4 * re);

    rhs_ch4_ice.x[i][j][k] =
        - transport_ch4_ice
        + diffusion_ch4_ice / (sc_ch4 * re);

    // Stokes terminal velocity for NH4SH crystals falling in the -r direction.
    // v_stokes [m/s] = (2/9) * r_p² * (rho_crystal - rho_mix) * g / mue_mix
    // Divided by u_0 to get the non-dimensional sedimentation velocity;
    // positive sign because downward settling ≡ negative radial velocity,
    // giving +v_sed * dq/dr in the concentration equation.
    const double v_stokes_nh4sh =
        (2.0 / 9.0) * r_p_nh4sh * r_p_nh4sh
        * (rg_nh4sh - rho_mix.x[i][j][k]) * g / mue_mix;

    rhs_nh4sh.x[i][j][k] =
        - transport_nh4sh
        + fluxlim_nh4sh.x[i][j][k]
        + diffusion_nh4sh / (sc_nh4sh * re)
        + chemical_reaction * massflux_nh4sh.x[i][j][k]
        + (v_stokes_nh4sh / u_0) * dnh4shdr;

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr_term;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe_term;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi_term;
}
/*
*
*/
