/*
 * Atmosphere General Circulation Modell(ATJUP) applied to turbulent flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with additional transport equations for the condensable species and for the
 * turbulent kinetic energy k* and its dissipation dis* (epsilon* or omega*)
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
 *
 * Renamed from RHS_Jup.cpp when k* and dis* became prognostic here, mirroring
 * ATOM_Precipitation/atmosphere/RHS_Atm_Turb.cpp: the closure in TurbulenceJup.h now only
 * supplies nue*, prod and the wall/ABL conditioning, while the two turbulence transport
 * equations themselves are assembled below and integrated by RungeKutta_Jup_Turb.cpp.
*/

#include "cJupiterModel.h"
#include "TurbulenceJup.h"   // shares nue_max_phys() with the closure

#include <cstdlib>   // getenv/atof for the radiative-coupling knob
#include <string>

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
    double dtkedr, ddisdr;

    double dudthe, dvdthe, dwdthe, dtdthe, dpdthe;
    double dh2odthe, dh2ocdthe, dh2oidthe;
    double dh2sdthe;
    double dnh3dthe, dnh3cdthe, dnh3idthe, dnh4shdthe;
    double dch4dthe, dch4cdthe, dch4idthe;
    double dtkedthe, ddisdthe;

    double dudphi, dvdphi, dwdphi, dtdphi, dpdphi;
    double dh2odphi, dh2ocdphi, dh2oidphi;
    double dh2sdphi;
    double dnh3dphi, dnh3cdphi, dnh3idphi, dnh4shdphi;
    double dch4dphi, dch4cdphi, dch4idphi;
    double dtkedphi, ddisdphi;

    // ---- Second-order derivative storage ----
    double d2udr2, d2vdr2, d2wdr2, d2tdr2;
    double d2h2odr2, d2h2ocdr2, d2h2oidr2;
    double d2h2sdr2;
    double d2nh3dr2, d2nh3cdr2, d2nh3idr2, d2nh4shdr2;
    double d2ch4dr2, d2ch4cdr2, d2ch4idr2;
    double d2tkedr2, d2disdr2;

    double d2udthe2, d2vdthe2, d2wdthe2, d2tdthe2;
    double d2h2odthe2, d2h2ocdthe2, d2h2oidthe2;
    double d2h2sdthe2;
    double d2nh3dthe2, d2nh3cdthe2, d2nh3idthe2, d2nh4shdthe2;
    double d2ch4dthe2, d2ch4cdthe2, d2ch4idthe2;
    double d2tkedthe2, d2disdthe2;

    double d2udphi2, d2vdphi2, d2wdphi2, d2tdphi2;
    double d2h2odphi2, d2h2ocdphi2, d2h2oidphi2;
    double d2h2sdphi2;
    double d2nh3dphi2, d2nh3cdphi2, d2nh3idphi2, d2nh4shdphi2;
    double d2ch4dphi2, d2ch4cdphi2, d2ch4idphi2;
    double d2tkedphi2, d2disdphi2;


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
    COMPUTE_DR(tke,       dtkedr,  d2tkedr2)
    COMPUTE_DR(dis,       ddisdr,  d2disdr2)
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
    COMPUTE_DTHE(tke,       dtkedthe,  d2tkedthe2)
    COMPUTE_DTHE(dis,       ddisdthe,  d2disdthe2)
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
    COMPUTE_DPHI(tke,       dtkedphi,  d2tkedphi2)
    COMPUTE_DPHI(dis,       ddisdphi,  d2disdphi2)
    dpdphi = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dphi;
    #undef COMPUTE_DPHI


    // ===== Neumann condition for k*/dis* at SeaMount faces =====
    // The centered stencils above read the zeroed obstacle cells (zero_land_cells in
    // TurbulenceJup.h sets tke = dis = 0 inside the SeaMount), which puts a large negative
    // Laplacian on the fluid cell touching the obstacle and drags k* — and hence nue* — to
    // zero exactly in the shear layer where the wake turbulence is generated. Replace the
    // solid neighbour by the local value so the diffusive flux through a solid face is zero.
    // Same treatment as ATOM's land-face branches in RHS_Atm_Turb.cpp, and identical to what
    // TurbulenceJup::compute_sources() already does when it forms its own gradients.
    {
        const bool solid_im1 = (SeaMount.x[i-1][j][k] == 1.0);
        const bool solid_ip1 = (SeaMount.x[i+1][j][k] == 1.0);
        const bool solid_jm1 = (SeaMount.x[i][j-1][k] == 1.0);
        const bool solid_jp1 = (SeaMount.x[i][j+1][k] == 1.0);
        const bool solid_km1 = (SeaMount.x[i][j][k-1] == 1.0);
        const bool solid_kp1 = (SeaMount.x[i][j][k+1] == 1.0);

        if(solid_im1 || solid_ip1){
            const double tke_im1 = solid_im1 ? tke.x[i][j][k] : tke.x[i-1][j][k];
            const double tke_ip1 = solid_ip1 ? tke.x[i][j][k] : tke.x[i+1][j][k];
            const double dis_im1 = solid_im1 ? dis.x[i][j][k] : dis.x[i-1][j][k];
            const double dis_ip1 = solid_ip1 ? dis.x[i][j][k] : dis.x[i+1][j][k];
            dtkedr   = (tke_ip1 - tke_im1) * inv_2dr * exp_rm;
            d2tkedr2 = (tke_ip1 - 2.0*tke.x[i][j][k] + tke_im1) * inv_dr2 * exp_2_rm;
            ddisdr   = (dis_ip1 - dis_im1) * inv_2dr * exp_rm;
            d2disdr2 = (dis_ip1 - 2.0*dis.x[i][j][k] + dis_im1) * inv_dr2 * exp_2_rm;
        }
        if(solid_jm1 || solid_jp1){
            const double tke_jm1 = solid_jm1 ? tke.x[i][j][k] : tke.x[i][j-1][k];
            const double tke_jp1 = solid_jp1 ? tke.x[i][j][k] : tke.x[i][j+1][k];
            const double dis_jm1 = solid_jm1 ? dis.x[i][j][k] : dis.x[i][j-1][k];
            const double dis_jp1 = solid_jp1 ? dis.x[i][j][k] : dis.x[i][j+1][k];
            dtkedthe   = (tke_jp1 - tke_jm1) * inv_2dthe;
            d2tkedthe2 = (tke_jp1 - 2.0*tke.x[i][j][k] + tke_jm1) * inv_dthe2;
            ddisdthe   = (dis_jp1 - dis_jm1) * inv_2dthe;
            d2disdthe2 = (dis_jp1 - 2.0*dis.x[i][j][k] + dis_jm1) * inv_dthe2;
        }
        if(solid_km1 || solid_kp1){
            const double tke_km1 = solid_km1 ? tke.x[i][j][k] : tke.x[i][j][k-1];
            const double tke_kp1 = solid_kp1 ? tke.x[i][j][k] : tke.x[i][j][k+1];
            const double dis_km1 = solid_km1 ? dis.x[i][j][k] : dis.x[i][j][k-1];
            const double dis_kp1 = solid_kp1 ? dis.x[i][j][k] : dis.x[i][j][k+1];
            dtkedphi   = (tke_kp1 - tke_km1) * inv_2dphi;
            d2tkedphi2 = (tke_kp1 - 2.0*tke.x[i][j][k] + tke_km1) * inv_dphi2;
            ddisdphi   = (dis_kp1 - dis_km1) * inv_2dphi;
            d2disdphi2 = (dis_kp1 - 2.0*dis.x[i][j][k] + dis_km1) * inv_dphi2;
        }
    }


    // ===== Nondimensionalisation of the body forces =====
    //
    // rhs_u is an acceleration in units of u_0^2/L, L = L_atm*1e3 m. Any term written from
    // physical constants therefore needs the factor that carries it into those units, and the
    // three body forces below each need a DIFFERENT one, because each is built from a different
    // combination of dimensional quantities:
    //
    //   Coriolis     2*Omega*u_phys       -> multiply by L/u_0     = 1400
    //   centrifugal  Omega^2*r_phys       -> multiply by L^2/u_0^2 = 1.96e6   (r_phys = rm*L)
    //   buoyancy     see the note below   -> multiply by 1e5*L/u_0^2 = 1.4e6
    //
    // Without them the terms are not "small", they are in the wrong unit system, and the model
    // has never felt any of them: Coriolis measures ~1e-4 against transport of order 1.
    //
    // These are computed from the model's own u_0 and L_atm rather than written as numbers, so
    // they follow the configuration instead of silently going stale when it changes.
    //
    // WHY ONE SWITCH AND NOT THREE KNOBS TO TASTE. The earlier attempt raised the buoyancy alone
    // to 1.4e6 and the radial velocity ran to 1535 m/s in 150 iterations; the conclusion drawn
    // was that the scale "cannot be switched on". That was the wrong conclusion from a right
    // measurement. Buoyancy at full strength with Coriolis still 1400x too weak is not a more
    // physical model, it is a NON-ROTATING one being convected: nothing in it can turn a vertical
    // plume into a balanced flow. On this planet the two belong to one balance and have to arrive
    // together. ATJUP_NONDIM=1 turns on all three; the individual ATJUP_ND_* switches exist for
    // attribution, not for production runs.
    //
    // The pressure gradient is deliberately NOT in this list. It needs no factor at all — see
    // p_dyn_to_bar() in cJupiterModel.h for why, and leave ATJUP_PGRAD_SCALE at 1.0.
    static const int nd_all = [](){ const char* e = getenv("ATJUP_NONDIM"); return e ? atoi(e) : 0; }();
    static const int nd_cor_on = [](){
        const char* e = getenv("ATJUP_ND_COR");  return e ? atoi(e) : -1; }();
    static const int nd_cent_on = [](){
        const char* e = getenv("ATJUP_ND_CENT"); return e ? atoi(e) : -1; }();
    static const int nd_buoy_on = [](){
        const char* e = getenv("ATJUP_ND_BUOY"); return e ? atoi(e) : -1; }();

    const double L_m = L_atm * 1.0e3;                       // shell thickness in metres
    const double nd_cor  = ((nd_cor_on  >= 0 ? nd_cor_on  : nd_all) != 0) ? L_m / u_0            : 1.0;
    const double nd_cent = ((nd_cent_on >= 0 ? nd_cent_on : nd_all) != 0) ? L_m * L_m / (u_0*u_0) : 1.0;
    const double nd_buoy = ((nd_buoy_on >= 0 ? nd_buoy_on : nd_all) != 0) ? 1.0e5 * L_m / (u_0*u_0) : 1.0;

    // ===== Coriolis and centrifugal forces =====
    double Coriolis_rad  = nd_cor * -2.0 * omega * sinthe * w_ijk;
    double Coriolis_the  = nd_cor * +2.0 * omega * costhe * w_ijk;
    double Coriolis_phi  = nd_cor * +2.0 * omega * (-costhe * v_ijk + sinthe * u_ijk);

    // Centrifugal acceleration = Omega^2 * s * s_hat, with s = r*sin(theta) the distance from
    // the rotation axis and s_hat = sin(theta)*e_r + cos(theta)*e_theta the unit vector pointing
    // AWAY from it. So
    //     a_r     = +Omega^2 * r * sin^2(theta)
    //     a_theta = +Omega^2 * r * sin(theta) * cos(theta)
    // Both were wrong. The code had Omega^2*r and Omega^2*r*|sin(theta)| entering with a MINUS,
    // so the force pointed toward the axis instead of away from it, the radial part had no
    // sin^2 at all (full strength at the poles, where it must vanish), and the meridional part
    // had |sin| in place of sin*cos, which is neither the right magnitude nor equator-directed.
    //
    // sinthe here is CLAMPED to >= 0.55 for the 1/sin^2 metric divisions; that floor has no
    // business in a body force, so the true sine is recovered from the cosine. theta runs 0..pi,
    // so sin(theta) >= 0 and the positive root is the right one. costhe now carries its proper
    // hemispheric sign (see cJupiterModel.h::costhe_abs), which is what makes a_theta point
    // toward the equator in BOTH hemispheres rather than southward everywhere.
    const double sinthe_true = sqrt(std::max(0.0, 1.0 - costhe * costhe));
    //
    // A WARNING about switching this one on, which no factor can fix. At full strength the
    // radial part is Omega^2*R = 2.17 m/s2 at the equator, 8.4 % of g, and the meridional part
    // peaks near 1 m/s2. On the real planet nothing has to balance those: they are absorbed into
    // the geopotential, and the answer is Jupiter's oblateness — the equator sits 4600 km further
    // from the centre than the poles, and the surfaces of constant effective gravity ARE that
    // shape. This model has a spherical grid, a spherical lower boundary and a constant radial g,
    // so the meridional part has nothing to work against and would drive a permanent, entirely
    // spurious pole-to-equator acceleration. The physical way to carry it is to fold it into an
    // effective gravity and never write it as a force at all. Measured below.
    double centrifugal_rad = nd_cent * omega * omega * rm * sinthe_true * sinthe_true;
    double centrifugal_the = nd_cent * omega * omega * rm * sinthe_true * costhe;

    double coeff_energy_p = u_0 * u_0 / (cp_mix * t_ref);


    // ===== Transport terms (advection) =====
    double v_invrm = v_ijk * inv_rm;
    double w_invrs = w_ijk * inv_rmsinthe;

    // ===== Compression work, and the adiabatic lapse rate it was missing =====
    //
    // coeff_energy_p * v.grad(p) is the correct nondimensional form of (1/(rho*cp)) Dp/Dt: with
    // p in the same nondimensional kinematic units p_dyn is stored in, u_0^2/(cp*t_ref) times
    // v.grad(p) reduces term for term to (L/(u_0*t_ref)) * (1/(rho*cp)) Dp/Dt. The form was
    // never the problem. The problem is WHICH pressure it was given.
    //
    // It was given p_dyn alone. But a parcel moving vertically does its work against the
    // HYDROSTATIC pressure, and that is the entire adiabatic lapse rate: with dp_stat/dz = -rho*g
    // the term reduces to dT/dt = -(g/cp)*w, which for Jupiter is 2.07 K per kilometre climbed.
    // Left out, the model has no adiabatic cooling at all, and therefore no static stability: a
    // parcel pushed up keeps the temperature it started with, arrives warmer than its new
    // surroundings, and is pushed up again. That is why the buoyancy could not be switched on.
    // It is not a scaling problem and no factor in rhs_u would have fixed it — with the buoyancy
    // off, nothing ever moved vertically for long enough to notice the term was absent.
    //
    // p_stat is in bar, so it is divided by p_dyn_to_bar() to enter as the same nondimensional
    // pressure. Only the radial derivative is taken: p_stat is a function of height alone here,
    // and its horizontal derivatives are zero by construction.
    //
    // It DEFAULTS TO ATJUP_NONDIM rather than to off, because a model that feels buoyancy without
    // it is not stably stratified and there is no configuration in which that is what you want.
    // ATJUP_ADIABATIC=0 forces it off even with the body forces on, which is how its contribution
    // was measured; =1 forces it on without them. With ATJUP_NONDIM unset it is off and the run is
    // bit-identical to before.
    static const bool adiabatic = [](){
        const char* e = getenv("ATJUP_ADIABATIC");
        if(e) return atoi(e) != 0;
        const char* n = getenv("ATJUP_NONDIM");
        return n && atoi(n) != 0; }();

    double dpstatdr = 0.0;
    if(adiabatic){
        const double to_nd = 1.0 / p_dyn_to_bar();      // bar -> nondimensional kinematic
        dpstatdr = (p_stat.x[i+1][j][k] - p_stat.x[i-1][j][k]) * inv_2dr * exp_rm * to_nd;
    }

    double pressure_t = coeff_energy_p
        * (u_ijk * (dpdr + dpstatdr) + v_invrm * dpdthe + w_invrs * dpdphi);

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

    // Advective form only (v.grad k, v.grad dis) — the standard k-eps / k-omega derivation
    // (Pope, Wilcox, Menter). ATOM once carried the conservative-form correction +k*div(v)
    // and dropped it again: the model is compressible, so in a persistently convergent region
    // (top of an updraft) that term grows k exponentially while its sink beta*.k.omega is only
    // linear in k. Do not reintroduce it here.
    double transport_tke       = u_ijk * dtkedr   + v_invrm * dtkedthe   + w_invrs * dtkedphi;
    double transport_dis       = u_ijk * ddisdr   + v_invrm * ddisdthe   + w_invrs * ddisdphi;


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

    double diffusion_tke = d2tkedr2 + dtkedr * two_inv_rm + d2tkedthe2 * inv_rm2
        + dtkedthe * costhe_inv_rm2sinthe + d2tkedphi2 * inv_rm2sinthe2;
    double diffusion_dis = d2disdr2 + ddisdr * two_inv_rm + d2disdthe2 * inv_rm2
        + ddisdthe * costhe_inv_rm2sinthe + d2disdphi2 * inv_rm2sinthe2;


    // ===== RHS assembly =====
    // The pressure gradient is in the same bind as the buoyancy: p_dyn is in bar, and rhs_u is
    // nondimensional in u_0^2/L_atm, so the consistent factor on (1/rho)*grad(p) is
    // 1e5/(r_mix*u_0^2) = 7.79. It matters only when the buoyancy is scaled, because those two
    // are the pair that must balance: raising the buoyancy alone by its own factor leaves
    // nothing able to oppose it, and the radial velocity runs away (measured: max|u| 38 -> 1535
    // m/s by iteration 150 with ATJUP_BUOY_SCALE=1.4e6 alone). Default 1.0 is bit-identical.
    static const double pgrad_scale = [](){
        const char* e = getenv("ATJUP_PGRAD_SCALE"); return e ? atof(e) : 1.0; }();

    double dpdr_term   = pgrad_scale * dpdr;
    double dpdthe_term = pgrad_scale * dpdthe * inv_rm;
    double dpdphi_term = pgrad_scale * dpdphi * inv_rmsinthe;

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
        // The density here used to be recomputed from a hard-coded R_H2He = 3600 J/(kg K),
        // while the model's own mixture constant, assembled from the composition in
        // ChemistryJup::ThermalPropertiesJup, is R_mix = 3104.84. Two different gas constants
        // for the same gas: the density came out 13.8 % too small and this heating 16 % too
        // large. rho_mix now carries exactly this quantity (computeMixtureDensity, refreshed at
        // the top of the physics block), so it is read rather than recomputed — one formula,
        // one constant. The ideal-gas fallback keeps the term alive if rho_mix is not yet filled.
        const double T_phys = t.x[i][j][k] * t_ref;             // [K]
        const double P_phys = p_stat.x[i][j][k] * 1.0e5;        // p_stat ~ bars -> [Pa]
        const double rho_f  = rho_mix.x[i][j][k];
        const double rho    = (rho_f > 0.0 && std::isfinite(rho_f))
                            ? rho_f
                            : ((T_phys > 1.0) ? P_phys / (R_mix * T_phys) : 0.0);  // [kg/m3]
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

    // ===== Turbulent (eddy) diffusion from the closure (opt-in) =====
    // nue is the DIMENSIONLESS eddy viscosity nue* = nue_phys/(u_0*L_atm), the same
    // normalisation as 1/re (since re = u_0*L_atm/nue_mol), so the two are directly additive:
    // the effective diffusivity is simply 1/re + nue*. Without this the closure computed an eddy
    // viscosity that nothing ever used. Scalars get nue*/Pr_t with a turbulent Prandtl (Schmidt)
    // number of 0.9, the standard value for shear-driven turbulence, in place of the laminar
    // 1/(sc*re). Gated by ATJUP_TURB_COUPLING (default 0 = off, bit-identical); nue is nonzero
    // only when ATJUP_TURB is enabled.
    static const double turb_coupling = [](){ const char* e = getenv("ATJUP_TURB_COUPLING"); return e ? atof(e) : 0.0; }();
    constexpr double Pr_t = 0.9;
    const double nue_t   = (turb_coupling != 0.0 && std::isfinite(nue.x[i][j][k]))
                         ? turb_coupling * std::max(0.0, nue.x[i][j][k]) : 0.0;
    const double nue_t_s = nue_t / Pr_t;      // scalar (heat / species) eddy diffusivity

    // Wall-adjacent eddy viscosity around the obstacle, filled once by computeWallViscosity()
    // (InitValues_Jup.cpp — the measurement that fixes its strength is documented there). It is
    // added to the MOMENTUM diffusivity only: the runaway is a momentum mode, and smearing the
    // condensable species across the obstacle face would be a chemistry change, not a numerical
    // one. Zero everywhere except within a few cells of the SeaMount.
    const double nue_wall = wall_nue.x[i][j][k];


    // ===== Turbulence closure: k* and dis* source terms and eddy diffusion coefficients =====
    // This is the ATOM arrangement (RHS_Atm_Turb.cpp): k* and dis* are two more transported
    // scalars of the Runge-Kutta system, so their production/destruction balance belongs here,
    // recomputed from the CURRENT tke/dis at every RK4 sub-stage. TurbulenceJup::run() only
    // executes on even iterations, so the stored tke_source/dis_source would otherwise be one
    // outer iteration stale, and constant across the four stages of the stage-splitting.
    //
    // The formulae are those of TurbulenceJup.h (which mirrors ATOM's TurbulenceAtm.h) — the
    // near-wall Chien damping of k-epsilon and the vorticity/strain correction of k-omega are
    // carried over too, so nothing the closure knows is dropped by moving the balance here.
    // What TurbulenceJup keeps is the eddy viscosity nue* itself, the friction velocity, the
    // wall/ABL conditioning and the land masking; nue* is read (not rewritten) below.
    //
    // NOTE on the diffusion coefficient: ATOM writes 1/re_turb + nue*/sigma, its re_turb being
    // the friction Reynolds number of the surface layer. ATJUP's RHS uses the model Reynolds
    // number re for the molecular background of every other equation, so k*/dis* use it too —
    // mixing the two would give the turbulence equations a molecular floor several orders of
    // magnitude away from that of the momentum equations they are coupled to.
    static const int turb_on = [](){ const char* e = getenv("ATJUP_TURB"); return e ? atoi(e) : 0; }();
    // 0 = k-epsilon (Chien 1982), 1 = k-omega (Wilcox 1988), 2 = k-omega SST (Menter 1994).
    // Same selection rule as TurbulenceJup: param.py's turb_model, overridable by ATJUP_TURB_MODEL,
    // anything unrecognised (including "none") falling back to SST.
    static const int turb_sel = [&](){
        const char* e = getenv("ATJUP_TURB_MODEL");
        const std::string s = e ? std::string(e) : turb_model;
        if(s == "k_epsilon") return 0;
        if(s == "k_omega")   return 1;
        return 2;
    }();

    double diffusion_tke_re = 0.0;
    double diffusion_dis_re = 0.0;
    double tke_src          = 0.0;
    double dis_src          = 0.0;

    if(turb_on != 0){
        constexpr double C_nue        = 0.028;      // C_mu, mirrors TurbulenceJup::C_nue
        constexpr double dis_min      = 1.0e-10;    // mirrors TurbulenceJup::dis_min
        constexpr double nue_gas_phys = 1.8e-5;     // kin. viscosity of the H2/He mix [m2/s]
        const double L_m       = L_atm * 1.0e3;     // L_atm is in km here, ATOM's is in metres
        const double nue_max   = TurbulenceJup::nue_max_phys() / (u_0 * L_m);
        const double nue_here  = std::isfinite(nue.x[i][j][k])
                               ? std::max(0.0, nue.x[i][j][k]) : 0.0;
        const double tke_s     = std::max(tke.x[i][j][k], dis_min);
        const double dis_s     = std::max(dis.x[i][j][k], dis_min);

        // Physical (per-metre, per-radian-arc) velocity gradients: the r-derivatives already
        // carry the stretching factor exp_rm from COMPUTE_DR, the angular ones still need the
        // 1/rm and 1/(rm sin(theta)) metric factors — the convention of compute_sources().
        const double dudr_s   = dudr;
        const double dvdr_s   = dvdr;
        const double dwdr_s   = dwdr;
        const double dudthe_s = dudthe * inv_rm;
        const double dvdthe_s = dvdthe * inv_rm;
        const double dwdthe_s = dwdthe * inv_rm;
        const double dudphi_s = dudphi * inv_rmsinthe;
        const double dvdphi_s = dvdphi * inv_rmsinthe;
        const double dwdphi_s = dwdphi * inv_rmsinthe;

        const double dtkedthe_s = dtkedthe * inv_rm;
        const double ddisdthe_s = ddisdthe * inv_rm;
        const double dtkedphi_s = dtkedphi * inv_rmsinthe;
        const double ddisdphi_s = ddisdphi * inv_rmsinthe;

        // grad(k*) . grad(dis*), the cross-diffusion contraction shared by k-omega and SST
        const double grad_dot = dtkedr * ddisdr
                              + dtkedthe_s * ddisdthe_s
                              + dtkedphi_s * ddisdphi_s;

        // ---- production tensor contraction P_k ----
        // cnue is rebuilt from the current tke/dis rather than read from nue.x, which lags a
        // sub-stage behind; identical to what compute_sources() does for the same reason.
        double cnue = (turb_sel == 0) ? C_nue * tke_s * tke_s / dis_s : tke_s / dis_s;
        cnue = std::min(cnue, nue_max);

        const double der = 0.66667 * (dudr_s + dvdthe_s + dwdphi_s);
        const double P_full = std::max(0.0,
              (cnue * (2.0 * dudr_s   - der) - 0.66667 * tke.x[i][j][k]) * dudr_s
            + (cnue * (dudthe_s + dvdr_s))                               * dudthe_s
            + (cnue * (dudphi_s + dwdr_s))                               * dudphi_s
            + (cnue * (2.0 * dvdthe_s - der) - 0.66667 * tke.x[i][j][k]) * dvdthe_s
            + (cnue * (dvdr_s + dudthe_s))                               * dvdr_s
            + (cnue * (dvdphi_s + dwdthe_s))                             * dvdphi_s
            + (cnue * (2.0 * dwdphi_s - der) - 0.66667 * tke.x[i][j][k]) * dwdphi_s
            + (cnue * (dwdr_s + dudphi_s))                               * dwdr_s
            + (cnue * (dwdthe_s + dvdphi_s))                             * dwdthe_s);
        prod.x[i][j][k] = P_full;

        // ---- wall distance above the local SeaMount surface, floored at one grid layer ----
        // TurbulenceJup uses the same floor: the blending/damping functions divide by y_star,
        // and the layer sitting on the surface has a wall distance of exactly zero.
        const double y_mount  = (double)get_layer_height(i_topography[j][k]) * 1.0e3;   // [m]
        const double dz_layer = std::max(1.0, layer_thickness_m(std::max(i - 1, 0)));
        const double y_phys   = std::max((double)get_layer_height(i) * 1.0e3 - y_mount, dz_layer);
        const double y_star   = y_phys / L_m;

        if(turb_sel == 0){                                   // k-epsilon, Chien 1982
            constexpr double sig_k   = 1.0;
            constexpr double sig_w   = 1.3;
            constexpr double C_eps_1 = 1.35;
            constexpr double C_eps_2 = 1.80;

            diffusion_tke_re = 1.0 / re + nue_here / sig_k;
            diffusion_dis_re = 1.0 / re + nue_here / sig_w;

            const double d_plus = y_phys * vel_star.y[j][k] / nue_gas_phys;
            const double Re_T   = tke_s * tke_s * u_0 * L_m / (dis_s * nue_gas_phys);
            const double f_2    = 1.0 - 0.4 / 1.8 * std::exp(-Re_T * Re_T / 36.0);

            // Chien near-wall damping, kept from TurbulenceJup::compute_k_epsilon
            const double y_star2 = y_star * y_star;
            const double L_k = -2.0 * tke.x[i][j][k] / y_star2;
            const double L_w = -2.0 * dis.x[i][j][k] / y_star2 * std::exp(-0.5 * d_plus);

            const double P_k = P_full;
            const double Y_k = dis.x[i][j][k];
            const double P_w = C_eps_1 * dis.x[i][j][k] / tke_s * P_k;
            const double Y_w = C_eps_2 * f_2 * dis.x[i][j][k] * dis.x[i][j][k] / tke_s;

            tke_src = P_k - Y_k + L_k / re_turb;
            dis_src = P_w - Y_w + L_w / re_turb;
        }
        else if(turb_sel == 1){                              // k-omega, Wilcox 1988/2006
            constexpr double sig_k    = 0.6;
            constexpr double sig_w    = 0.5;
            constexpr double bet_star = 0.09;
            constexpr double gam      = 0.52;
            constexpr double bet_0    = 0.0708;
            constexpr double C_lim    = 0.875;

            // Wilcox multiplies by sigma where k-epsilon and SST divide — the standard form
            // of each model, mirrored from ATOM.
            diffusion_tke_re = 1.0 / re + sig_k * nue_here;
            diffusion_dis_re = 1.0 / re + sig_w * nue_here;

            const double S11 = dudr_s, S22 = dvdthe_s, S33 = dwdphi_s;
            const double S12 = 0.5 * (dudthe_s + dvdr_s);
            const double S13 = 0.5 * (dudphi_s + dwdr_s);
            const double S23 = 0.5 * (dvdphi_s + dwdthe_s);
            const double S_mag = std::sqrt(2.0 * (S11*S11 + S22*S22 + S33*S33
                                                + 2.0*(S12*S12 + S13*S13 + S23*S23)));

            const double W12 = dudthe_s - dvdr_s;
            const double W13 = dudphi_s - dwdr_s;
            const double W23 = dvdphi_s - dwdthe_s;
            const double Omega = std::sqrt(W12*W12 + W13*W13 + W23*W23);

            // Wilcox (2006) vorticity-strain correction of the omega destruction coefficient
            const double chi_w  = std::fabs(Omega * Omega * S_mag
                                / std::pow(bet_star * dis_s, 3));
            const double f_bet  = (1.0 + 85.0 * chi_w) / (1.0 + 100.0 * chi_w);
            const double bet_wc = bet_0 * f_bet;

            const double sig_d = (grad_dot <= 0.0) ? 0.0 : 0.125;

            const double P_k = std::min(P_full, 20.0 * bet_star * tke_s * dis_s);
            const double Y_k = bet_star * tke_s * dis_s;
            const double P_w = gam * dis_s / tke_s * P_k;
            const double Y_w = bet_wc * dis_s * dis_s;
            const double D_w = sig_d / dis_s * grad_dot;

            tke_src = P_k - Y_k;
            dis_src = P_w - Y_w + D_w;

            // C_lim is the stress limiter of the 2006 revision; it acts on nue*, which
            // TurbulenceJup owns, so it is only referenced here to keep the constant set
            // complete and identical between the two files.
            (void)C_lim;
        }
        else {                                               // k-omega SST, Menter 1994
            constexpr double bet_star = 0.09;   // beta* destruction coefficient, NOT C_mu
            constexpr double sig_k1   = 1.176;
            constexpr double sig_k2   = 1.0;
            constexpr double sig_w1   = 2.0;
            constexpr double sig_w2   = 1.168;
            constexpr double bet1     = 0.0333;
            constexpr double bet2     = 0.0368;
            constexpr double gam1     = 0.413;
            constexpr double gam2     = 0.2;

            const double nue_air_nd = nue_gas_phys / (u_0 * L_m);
            const double CD_kw = std::max(2.0 * sig_w2 / dis_s * grad_dot, 1.0e-20);

            // Menter's arg1: the molecular viscosity, not the turbulent one, enters the
            // 500 nue/(y^2 omega) branch.
            const double arg1 = std::min(
                std::max(std::sqrt(std::max(tke.x[i][j][k], 0.0)) / (bet_star * dis_s * y_star),
                         500.0 * nue_air_nd / (y_star * y_star * dis_s)),
                4.0 * sig_w2 * tke_s / (CD_kw * y_star * y_star));
            const double F1 = std::tanh(std::pow(arg1, 4));

            // blend(inner, outer, F1) = F1*inner + (1-F1)*outer, as in TurbulenceJup
            const double sig_k = F1 * sig_k1 + (1.0 - F1) * sig_k2;
            const double sig_w = F1 * sig_w1 + (1.0 - F1) * sig_w2;
            diffusion_tke_re = 1.0 / re + nue_here / sig_k;
            diffusion_dis_re = 1.0 / re + nue_here / sig_w;

            const double P_k = std::min(P_full, 20.0 * bet_star * tke_s * dis_s);
            const double Y_k = bet_star * tke_s * dis_s;
            const double P_w = (F1 * gam1 + (1.0 - F1) * gam2) * P_k * dis_s / tke_s;
            const double Y_w = (F1 * bet1 + (1.0 - F1) * bet2) * dis_s * dis_s;
            const double D_w = 2.0 * (1.0 - F1) * sig_w2 / dis_s * grad_dot;

            tke_src = P_k - Y_k;
            dis_src = P_w - Y_w + D_w;
        }

        // Same source caps as compute_sources(): they bound the P-Y balance and stop the
        // 1/sin^2(theta) amplification of the cross-diffusion term blowing up at the poles.
        const double tke_src_max = 20.0 * tke_s * dis_s;
        const double dis_src_max = 20.0 * dis_s * dis_s;
        if(!std::isfinite(tke_src)) tke_src = 0.0;
        if(!std::isfinite(dis_src)) dis_src = 0.0;
        tke_src = std::max(-tke_src_max, std::min(tke_src_max, tke_src));
        dis_src = std::max(-dis_src_max, std::min(dis_src_max, dis_src));

        // Publish for ParaView / printMinMax, exactly as ATOM's RHS does.
        tke_source.x[i][j][k] = tke_src;
        dis_source.x[i][j][k] = dis_src;
    }

    rhs_t.x[i][j][k] =
        + pressure_t
        - transport_t
        + diffusion_t * (1.0 / (re * pr) + nue_t_s)
        + radiation_t
        + precip_t;

    // ===== Buoyancy: a Boussinesq anomaly, with the sign and the scale it needs =====
    //
    // rho = p/(R*T), and buoy_ref_level[i] (computeBuoyancyRefLevel, once per RK4 step) is the
    // area-weighted horizontal mean of the same expression at this level, so the difference is
    // (g/r_mix)*(rho - rho_bar): zero mean at every height, only horizontal density contrasts
    // drive vertical motion, and hydrostatic balance is left to carry the mean.
    //
    // TWO THINGS WERE WRONG WITH IT, and they had to be fixed together.
    //
    // (1) THE SIGN. The anomaly entered rhs_u with a PLUS, so a parcel DENSER than its level
    // mean was accelerated UPWARD and a warm, light one pushed down — convection upside down.
    // The Archimedes force is a = -g*(rho - rho_bar)/rho_ref. The error came from translating
    // ATOM's form, which is +(t - t_ref_level[i]): a plus is right for a TEMPERATURE anomaly,
    // because warm rises, and wrong for a DENSITY anomaly, because heavy sinks. The sign flip
    // that the change of variable requires was dropped.
    //
    // (2) THE SCALE. p_stat is in bar while r_mix*R_mix*T yields pascals, so the expression is
    // 1e-5 of a physical acceleration; and rhs_u is nondimensional in units of u_0^2/L_atm =
    // 0.0714 m/s2, so a physical acceleration still needs L_atm/u_0^2 = 14. The dimensionally
    // consistent factor is therefore 1e5 * L_atm/u_0^2 = 1.4e6. Forces() carries only the 1e5,
    // and correctly so — BuoyancyForce there is a diagnostic force DENSITY in N/m3, matching
    // CoriolisForce, a different unit system from this equation.
    //
    // Fixing (2) without (1) would have been the worst of both: the term is currently ~2e-4
    // against O(1) transport, i.e. switched off, so the wrong sign costs nothing today. Scale
    // it up by 1.4e6 with the sign still inverted and the model would convect upside down at
    // full strength.
    //
    // The sign is corrected unconditionally, because it is a defect. The scale is a knob,
    // because g*L_atm/u_0^2 = 363 against O(1) transport is a violent change to a model that
    // has never felt buoyancy, and ATOM's experience with the same correction (a 336x factor,
    // see cAtmosphereModel.h) was that switching it on cold blows the CFL limit and needs a
    // ramp of a few hundred iterations:
    //   ATJUP_NONDIM=1                 the dimensionally consistent value, together with the
    //                                  Coriolis and centrifugal factors it belongs with
    //   ATJUP_BUOY_SCALE=<x>           a MULTIPLIER on top of that, for finding the usable
    //                                  range; 1.0 = the physical value. (It used to carry the
    //                                  1.4e6 itself — do not write that number here any more,
    //                                  ATJUP_NONDIM supplies it and the two would multiply.)
    //   ATJUP_BUOY_RAMP_ITERS=<n>      ramp linearly from 0 to full over the first n iterations
    //
    // NOT addressed, and worth knowing before the knob is turned up: the anomaly is divided by
    // the CONSTANT r_mix, not by the level mean density. Proper Boussinesq divides by the local
    // reference rho_bar(i), which would weight the anomaly aloft up to 200x more strongly than
    // this does. rho_mix and buoy_ref_level are both available if that is wanted.
    //
    // ---- Why the density here is built from p_stat alone (ATJUP_BUOY_PDYN=1 puts p_dyn back) ----
    //
    // The buoyancy used to read rho = (p_stat + p_dyn)/(R*T), which closes a loop that has no
    // physics in it: buoyancy drives a divergence, the divergence sets p_dyn, p_dyn changes the
    // density, and the density feeds the buoyancy again. In a Boussinesq or anelastic system the
    // density anomaly is a THERMODYNAMIC quantity — it comes from temperature and composition
    // against a hydrostatic reference pressure — while p_dyn is a Lagrange multiplier enforcing
    // the velocity constraint. The two are not the same kind of object and the second does not
    // belong in the equation of state.
    //
    // The loop was invisible as long as p_dyn was a local smear of the divergence: with one
    // Gauss-Seidel sweep per step p_dyn stays near 0.002 bar against a p_stat of several bar. It
    // becomes fatal the moment the elliptic problem is actually solved. Measured, 99 iterations,
    // buoyancy at full scale: 1 sweep max|u| = 277 m/s, 50 sweeps 127, and at 300 sweeps the run
    // reaches 69775 m/s by iteration 33 and collapses, with p_dyn at +-60 bar. More convergence
    // made it worse, which is the signature of a feedback rather than of a discretisation error.
    static const double buoy_scale = [](){
        const char* e = getenv("ATJUP_BUOY_SCALE"); return e ? atof(e) : 1.0; }();
    static const int buoy_ramp_iters = [](){
        const char* e = getenv("ATJUP_BUOY_RAMP_ITERS"); return e ? atoi(e) : 0; }();
    const double buoy_ramp = (buoy_ramp_iters > 0)
        ? std::min(1.0, (double)iter_n / (double)buoy_ramp_iters) : 1.0;

    // The 1/t is guarded because a cell without a positive temperature has no density and
    // hence no buoyancy: t = 0 used to give +-inf here, RK4 turned that into an infinite
    // velocity, and the advection stencil of every neighbour then carried inf - inf = NaN
    // outwards, roughly doubling the affected volume each iteration.
    const double buoyancy_u = (t.x[i][j][k] > 0.0)
        ? -nd_buoy * buoy_scale * buoy_ramp * buoyancy
          * (g * buoy_pressure(i, j, k)
               / (r_mix * R_mix * t.x[i][j][k] * t_ref)
             - buoy_ref_level[i])
        : 0.0;

    rhs_u.x[i][j][k] =
        - dpdr_term
        + buoyancy_u
        - transport_u
        + diffusion_u * (1.0 / re + nue_t + nue_wall)
        - Coriolis    * Coriolis_rad
        + centrifugal * centrifugal_rad;


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
        + diffusion_v * (1.0 / re + nue_t + nue_wall)
        - Coriolis    * Coriolis_the
        + centrifugal * centrifugal_the;

    rhs_w.x[i][j][k] =
        - dpdphi_term
        - transport_w
        + diffusion_w * (1.0 / re + nue_t + nue_wall)
        - Coriolis    * Coriolis_phi;

    // ---- Single-cell momentum budget probe (ATJUP_PROBE="i,j,k", off by default) ----
    // The secular growth of max|w| is anchored to ONE fluid cell beside the staircase flank of
    // the SeaMount cone, so the question "which term feeds it" is answerable by printing the
    // terms of that one cell. Fires on all four RK4 stages, which also shows whether a term is
    // steady through the step or oscillating between stages. Velocities are in m/s, the terms
    // are left nondimensional: multiplied by dt they are the increment per stage.
    {
        static const bool probe_on = getenv("ATJUP_PROBE") != nullptr;
        if(probe_on){
            static int pi = -1, pj = -1, pk = -1;
            static bool parsed = false;
            if(!parsed){ sscanf(getenv("ATJUP_PROBE"), "%d,%d,%d", &pi, &pj, &pk); parsed = true; }
            if(i == pi && j == pj && k == pk){
                // diff is the FULL diffusion actually used, wall viscosity included, so the
                // four terms sum to rhs_w and the budget can be checked by eye.
                printf("      PROBE %4d  w=%9.3f u=%9.3f v=%9.3f | rhs_w=%11.4f"
                       " dpdphi=%11.4f transp=%11.4f diff=%11.4f cor=%11.4f"
                       " | nue_t=%.3e nue_wall=%.3e\n",
                       iter_n, w_ijk * u_0, u_ijk * u_0, v_ijk * u_0, rhs_w.x[i][j][k],
                       -dpdphi_term, -transport_w,
                       diffusion_w * (1.0 / re + nue_t + nue_wall),
                       -Coriolis * Coriolis_phi, nue_t, nue_wall);
                // Advection broken into its four pieces: the three directional derivatives and
                // the spherical metric group. Which one carries the +4 tells the difference
                // between "the flow really accelerates round the flank" (r/theta/phi advection)
                // and "the curvature terms are unbalanced at the wall" (metric).
                // Temperature budget of the same cell. dT/dt is nondimensional; multiplied by
                // dt*t_ref it is the temperature change per stage in kelvin, which is the form
                // to compare against the adiabatic lapse rate when ATJUP_ADIABATIC is on.
                printf("      PROBET %4d  T=%8.3fK rhs_t=%11.4e | compr=%11.4e"
                       " (dyn=%11.4e stat=%11.4e) transp=%11.4e diff=%11.4e rad=%11.4e\n",
                       iter_n, t.x[i][j][k] * t_ref, rhs_t.x[i][j][k], pressure_t,
                       coeff_energy_p * (u_ijk * dpdr + v_invrm * dpdthe + w_invrs * dpdphi),
                       coeff_energy_p * u_ijk * dpstatdr, -transport_t,
                       diffusion_t * (1.0 / (re * pr) + nue_t_s), radiation_t);
                printf("      PROBEADV %4d  u*dwdr=%11.4f  v/r*dwdthe=%11.4f"
                       "  w/(r sin)*dwdphi=%11.4f  metric=%11.4f | dwdr=%10.3f dwdthe=%10.3f"
                       " dwdphi=%10.3f\n",
                       iter_n, u_ijk * dwdr, v_invrm * dwdthe, w_invrs * dwdphi,
                       (w_ijk * u_ijk + v_ijk * w_ijk * cotanthe) * inv_rm,
                       dwdr, dwdthe, dwdphi);
            }
        }
    }

    rhs_h2o.x[i][j][k] =
        - transport_h2o
        + diffusion_h2o * (1.0 / (sc_h2o * re) + nue_t_s);

    rhs_h2o_cloud.x[i][j][k] =
        - transport_h2o_cloud
        + diffusion_h2o_cloud * (1.0 / (sc_h2o * re) + nue_t_s);

    rhs_h2o_ice.x[i][j][k] =
        - transport_h2o_ice
        + diffusion_h2o_ice * (1.0 / (sc_h2o * re) + nue_t_s);

    rhs_h2s.x[i][j][k] =
        - transport_h2s
        + diffusion_h2s * (1.0 / (sc_h2s * re) + nue_t_s)
        + chemical_reaction * massflux_h2s.x[i][j][k];

    rhs_nh3.x[i][j][k] =
        - transport_nh3
        + diffusion_nh3 * (1.0 / (sc_nh3 * re) + nue_t_s)
        + chemical_reaction * massflux_nh3.x[i][j][k];

    rhs_nh3_cloud.x[i][j][k] =
        - transport_nh3_cloud
        + diffusion_nh3_cloud * (1.0 / (sc_nh3 * re) + nue_t_s);

    rhs_nh3_ice.x[i][j][k] =
        - transport_nh3_ice
        + diffusion_nh3_ice * (1.0 / (sc_nh3 * re) + nue_t_s);

    rhs_ch4.x[i][j][k] =
        - transport_ch4
        + diffusion_ch4 * (1.0 / (sc_ch4 * re) + nue_t_s);

    rhs_ch4_cloud.x[i][j][k] =
        - transport_ch4_cloud
        + diffusion_ch4_cloud * (1.0 / (sc_ch4 * re) + nue_t_s);

    rhs_ch4_ice.x[i][j][k] =
        - transport_ch4_ice
        + diffusion_ch4_ice * (1.0 / (sc_ch4 * re) + nue_t_s);

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
        + diffusion_nh4sh * (1.0 / (sc_nh4sh * re) + nue_t_s)
        + chemical_reaction * massflux_nh4sh.x[i][j][k]
        + (v_stokes_nh4sh / u_0) * dnh4shdr;

    // ===== Turbulence transport equations =====
    // dk*/dt   = -v.grad k*   + div((1/re + nue*/sigma_k) grad k*)   + (P_k - Y_k)
    // ddis*/dt = -v.grad dis* + div((1/re + nue*/sigma_w) grad dis*) + (P_w - Y_w + D_w)
    // With the closure off both are identically zero, so RungeKuttaJup leaves k*/dis* at their
    // initial values and the run stays bit-identical to the pre-turbulence model.
    if(turb_on != 0){
        rhs_tke.x[i][j][k] =
            - transport_tke
            + diffusion_tke * diffusion_tke_re
            + tke_src;

        rhs_dis.x[i][j][k] =
            - transport_dis
            + diffusion_dis * diffusion_dis_re
            + dis_src;
    } else {
        rhs_tke.x[i][j][k] = 0.0;
        rhs_dis.x[i][j][k] = 0.0;
    }

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr_term;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe_term;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi_term;
}
/*
*
*/
