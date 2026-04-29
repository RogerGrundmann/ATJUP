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

    double dudthe, dvdthe, dwdthe, dtdthe, dpdthe;
    double dh2odthe, dh2ocdthe, dh2oidthe;
    double dh2sdthe;
    double dnh3dthe, dnh3cdthe, dnh3idthe, dnh4shdthe;

    double dudphi, dvdphi, dwdphi, dtdphi, dpdphi;
    double dh2odphi, dh2ocdphi, dh2oidphi;
    double dh2sdphi;
    double dnh3dphi, dnh3cdphi, dnh3idphi, dnh4shdphi;

    // ---- Second-order derivative storage ----
    double d2udr2, d2vdr2, d2wdr2, d2tdr2;
    double d2h2odr2, d2h2ocdr2, d2h2oidr2;
    double d2h2sdr2;
    double d2nh3dr2, d2nh3cdr2, d2nh3idr2, d2nh4shdr2;

    double d2udthe2, d2vdthe2, d2wdthe2, d2tdthe2;
    double d2h2odthe2, d2h2ocdthe2, d2h2oidthe2;
    double d2h2sdthe2;
    double d2nh3dthe2, d2nh3cdthe2, d2nh3idthe2, d2nh4shdthe2;

    double d2udphi2, d2vdphi2, d2wdphi2, d2tdphi2;
    double d2h2odphi2, d2h2ocdphi2, d2h2oidphi2;
    double d2h2sdphi2;
    double d2nh3dphi2, d2nh3cdphi2, d2nh3idphi2, d2nh4shdphi2;


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

    double diffusion_nh4sh = d2nh4shdr2 + dnh4shdr * two_inv_rm + d2nh4shdthe2 * inv_rm2
        + dnh4shdthe * costhe_inv_rm2sinthe + d2nh4shdphi2 * inv_rm2sinthe2;


    // ===== RHS assembly =====
    double dpdr_term   = dpdr;
    double dpdthe_term = dpdthe * inv_rm;
    double dpdphi_term = dpdphi * inv_rmsinthe;

    rhs_t.x[i][j][k] =
        + pressure_t
        - transport_t
        + diffusion_t / (re * pr);

    rhs_u.x[i][j][k] =
        - dpdr_term
        + buoyancy * g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])
                      / (r_mix * R_mix * t.x[i][j][k] * t_ref)
        - transport_u
        + diffusion_u / re
        - Coriolis    * Coriolis_rad
        - centrifugal * centrifugal_rad;

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

    rhs_nh4sh.x[i][j][k] =
        - transport_nh4sh
        + diffusion_nh4sh / (sc_nh4sh * re)
        + chemical_reaction * massflux_nh4sh.x[i][j][k];

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr_term;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe_term;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi_term;
}
/*
*
*/
