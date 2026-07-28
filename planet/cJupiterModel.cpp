/*
 * Jupiter Atmosphere Circulation Model (JACM) applied to laminar flow
 * program for the computation of jupiter-atmospherical circulating flows in a spherical shell
 * modeling of the atmosphere with gases, their cloud and ice formation: H2O, H2S, NH3 and NH4SH
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)
*/

#include "cJupiterModel.h"

#include "ChemistryJup.h"
#include "PressureSolverJup.h"
#include "SaturationAdjustmentJup.h"
#include "BC_Jup.h"
#include "VelocityInitializerJup.h"
#include "RadiationJup.h"
#include "PrecipitationJup.h"
#include "TurbulenceJup.h"
#include "ConvectiveAdjustmentJup.h"
#include "ThermalWindJup.h"

#include <cstdlib>   // getenv/atof/atoi for the Shapiro velocity-filter knobs

using namespace std;
using namespace tinyxml2;
using namespace JupiterUtils;

// ----------------------------------------------------------------------------
// Velocity Shapiro-filter knobs, read once from the environment.
//   ATJUP_VEL_SHAPIRO_ORDER   2 = 1-2-1 (default, bit-identical) | 4 = shear-preserving
//   ATJUP_SHAPIRO_STRENGTH    filter strength (default 1.0)
//   ATJUP_VEL_SHAPIRO_INLOOP  per-iteration passes applied to u,v,w after RK4 (default 0 = off)
// Defaults reproduce the previous behaviour bit-for-bit (order 2, strength 1, no in-loop pass).
static int    shapiro_vel_order()  { static const int    v = [](){ const char* e = getenv("ATJUP_VEL_SHAPIRO_ORDER");  return e ? atoi(e) : 2;   }(); return v; }
static double shapiro_strength()   { static const double v = [](){ const char* e = getenv("ATJUP_SHAPIRO_STRENGTH");   return e ? atof(e) : 1.0; }(); return v; }
static int    shapiro_vel_inloop() { static const int    v = [](){ const char* e = getenv("ATJUP_VEL_SHAPIRO_INLOOP"); return e ? atoi(e) : 0;   }(); return v; }

// Radiation scaffold knob (default off = bit-identical). When set, RadiationJup runs
// each even iteration and fills the diagnostic radiation/epsilon/Q_rad arrays; it does
// not yet feed the temperature equation.
static int    radiation_enabled()  { static const int    v = [](){ const char* e = getenv("ATJUP_RADIATION");           return e ? atoi(e) : 0;   }(); return v; }

// Precipitation microphysics knob (default off = bit-identical). When set, PrecipitationJup
// runs once per iteration and fills the rain/snow/graupel fluxes for H2O and NH3 plus the
// NH4SH settling flux, and the latent-heat field Q_precip. It DOES feed back: the condensate
// it converts is removed from the cloud/ice fields in place. That happens after the three
// SaturationAdjustmentJup calls, so the removal persists and the next iteration must draw on
// the vapour reservoir to rebuild cloud — which is what makes the rate self-limiting.
// Q_precip only reaches rhs_t if ATJUP_PRECIP_COUPLING is also set (see RHS_Jup_Turb.cpp).
// Turbulence knob (default off = bit-identical). Mirrors ATOM's TurbulenceAtm: all three models
// (k_epsilon | k_omega | k_omega_SST) are available, selected by cJupiterModel::turb_model or
// the ATJUP_TURB_MODEL environment variable. Fills tke/dis/nue/prod/tke_source/dis_source.
// With it on, k* and dis* are prognostic: RHS_Jup_Turb.cpp assembles rhs_tke/rhs_dis and
// RungeKutta_Jup_Turb.cpp integrates them, as in ATOM. nue* only reaches the momentum and
// scalar equations if ATJUP_TURB_COUPLING is also set.
static int    turb_enabled()       { static const int    v = [](){ const char* e = getenv("ATJUP_TURB");                 return e ? atoi(e) : 0;   }(); return v; }
static int    precip_enabled()     { static const int    v = [](){ const char* e = getenv("ATJUP_PRECIP");               return e ? atoi(e) : 0;   }(); return v; }

// Dry convective adjustment (ConvectiveAdjustmentJup). It DEFAULTS TO ATJUP_NONDIM rather than to
// off, because the two belong together: switching the buoyancy on without it gives a model that
// releases its superadiabatic layer and has no way to relieve it, and that run dies at iteration
// 233. ATJUP_CONV_ADJ=0 forces it off even with the body forces on, which is how the measurement
// of what it is worth was made; ATJUP_CONV_ADJ=1 forces it on without them, which shows what the
// stratification does on its own. With ATJUP_NONDIM unset the default is off and every run is
// bit-identical to before.
static int    conv_adj_enabled(){
    static const int v = [](){
        const char* e = getenv("ATJUP_CONV_ADJ");
        if(e) return atoi(e);
        const char* n = getenv("ATJUP_NONDIM");
        return n ? atoi(n) : 0;
    }();
    return v;
}

// Thermal-wind initialisation (ThermalWindJup), once, at the end of the initial state. Also
// defaults to ATJUP_NONDIM: with the body forces off there is no buoyancy for the wind to
// balance and the adjustment would be a change for nothing, while with them on an unbalanced
// initial state is the single largest thing wrong with the model. ATJUP_THERMAL_WIND=0/1 forces
// it either way.
static int    thermal_wind_enabled(){
    static const int v = [](){
        const char* e = getenv("ATJUP_THERMAL_WIND");
        if(e) return atoi(e);
        const char* n = getenv("ATJUP_NONDIM");
        return n ? atoi(n) : 0;
    }();
    return v;
}


cJupiterModel* cJupiterModel::m_model = NULL;

const double cJupiterModel::pi180 = 180.0/M_PI;      // pi180 = 57.3

const double cJupiterModel::the_degree = 1.0;         // compares to 1° step size laterally
const double cJupiterModel::phi_degree = 1.0;         // compares to 1° step size longitudinally

const double cJupiterModel::the0 = 0.0;             // North Pole
const double cJupiterModel::phi0 = 0.0;             // zero meridian in Greenwich

const double cJupiterModel::r0 = 1.0; // value much too small, Jupiter radius 72000km

const double cJupiterModel::dr = 0.025;    // 0.025 x 40 = 1.0 compares to 16 km : 40 = 400 m for 1 radial step
const double cJupiterModel::dthe = the_degree/pi180; 
const double cJupiterModel::dphi = phi_degree/pi180;


cJupiterModel::cJupiterModel():
    j_ellipse(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    has_welcome_msg_printed(false){
//    if(PythonStream::is_enable()){
//        backup = std::cout.rdbuf();
//        std::cout.rdbuf(&ps);
//    }
    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);
    // set default configuration
    SetDefaultConfig();
    m_model = this;
    rad.initArray_1D(im, 0); // radial coordinate direction
    the.initArray_1D(jm, 0); // lateral coordinate direction
    phi.initArray_1D(km, 0); // longitudinal coordinate direction
    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);
    init_layer_heights();
}

cJupiterModel::~cJupiterModel(){
    delete m_chem;
    m_chem = nullptr;

    delete m_pressure;
    m_pressure = nullptr;

    if(PythonStream::is_enable()){
        std::cout.rdbuf(backup);
    }
    m_model = NULL;
}

ChemistryJup& cJupiterModel::getChemistry(){
    if(!m_chem)
        m_chem = new ChemistryJup(*this);
    return *m_chem;
}

void cJupiterModel::ChemMassRateJup()     { getChemistry().ChemMassRateJup(); }
void cJupiterModel::DiffMassFluxJup()     { getChemistry().DiffMassFluxJup(); }
void cJupiterModel::ThermalPropertiesJup(){ getChemistry().ThermalPropertiesJup(); }

PressureSolverJup& cJupiterModel::getPressureSolver(){
    if(!m_pressure)
        m_pressure = new PressureSolverJup(*this);
    return *m_pressure;
}

// Does the radius the metric uses equal the radius Jupiter has?
//
// Ported from ATOM_Precipitation (aee862d/0c19e06), where the same ambiguity cost a working day.
// The point is not that the answer is unknown — it is written in a comment at r0's definition
// ("value much too small, Jupiter radius 72000km") and in several notes — but that a run gives no
// sign of which state it is in. RungeKuttaJup takes geo.rm = rad.z[i] and uses it as the PLANETARY
// radius in v/r d/dthe, 1/(r sinthe) d/dphi and every Laplacian, and rad.z runs 1.000 .. 2.000
// unless ATJUP_METRIC_RADIUS moves it. One rad.z unit is L_atm = 140 km (see
// metricShellLength_km), so the metric puts Jupiter's surface at 140 km from the centre instead
// of 69911, and every horizontal derivative is 499x too large.
//
// It WARNS and continues, because that is the model's default state rather than a regression, and
// a hard stop would prevent every unconverted run from starting. ATJUP_METRIC_STRICT=1 makes it
// fatal — for use once ATJUP_METRIC_RADIUS is the default, so the two cannot drift apart again.
void cJupiterModel::checkMetricConsistency() const {
    const double L_unit_km  = metricShellLength_km();      // km per rad.z unit
    const double r_metric_km = rad.z[0] * L_unit_km;       // radius the metric implies
    const double ratio = (r_metric_km > 0.0) ? R_planet_km / r_metric_km : 0.0;
    const bool consistent = std::fabs(ratio - 1.0) < 1.0e-3;

    printf("\n      ATJUP: metric check - core radius = %.1f km (rad.z[0] %.3f x %.1f km per unit),"
           " planet radius = %.1f km\n", r_metric_km, rad.z[0], L_unit_km, R_planet_km);

    if(consistent){
        printf("      ATJUP: metric check - consistent.\n");
        return;
    }

    printf("      ATJUP: metric check - MISMATCH by a factor of %.3f.\n", ratio);
    printf("            Horizontal metric terms (inv_rm, inv_rmsinthe in RungeKutta_Jup_Turb,\n"
           "            PressureSolverJup, TurbulenceJup) are that factor too large.\n");
    printf("            Set ATJUP_METRIC_RADIUS=%.0f to correct it.\n", R_planet_km);

    static const bool strict = [](){
        const char* e = getenv("ATJUP_METRIC_STRICT"); return e && atoi(e) != 0; }();
    if(strict){
        printf("      ATJUP: metric check - ATJUP_METRIC_STRICT is set, stopping.\n");
        std::exit(1);
    }
}

// Shapiro de-checkerboarding of the velocity fields. Order 4 preserves the resolved
// zonal-jet shear (∂w/∂θ) and the GRS wake far better than the 1-2-1 (order 2) while
// still annihilating the 2Δ grid mode; order 2 reproduces the previous behaviour.
void cJupiterModel::dampVelocities(){
    const double s = shapiro_strength();
    if(shapiro_vel_order() == 4){
        JupiterUtils::damp_wiggles_ho(u, &i_topography, true, true, true, s);
        JupiterUtils::damp_wiggles_ho(v, &i_topography, true, true, true, s);
        JupiterUtils::damp_wiggles_ho(w, &i_topography, true, true, true, s);
    }else{
        JupiterUtils::damp_wiggles(u, &i_topography, true, true, true, s);
        JupiterUtils::damp_wiggles(v, &i_topography, true, true, true, s);
        JupiterUtils::damp_wiggles(w, &i_topography, true, true, true, s);
    }
}


#include "cJupiterDefaults.cpp.inc"
/*
*
*/
void cJupiterModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    if(err){
        doc.PrintError();
        throw std::invalid_argument("   couldn't load config file inside cJupiterModel");
    }
    XMLElement *atjup = doc.FirstChildElement("atjup");
    if(!atjup){
        return;
    }
    XMLElement* elem_common = doc.FirstChildElement("atjup")->FirstChildElement("common");
    if(!elem_common){
        return;
    }
    XMLElement* elem_jupiter = doc.FirstChildElement("atjup")->FirstChildElement("jupiter");
    if(!elem_jupiter){
        return;
    }
#include "JupiterLoadConfig.cpp.inc"
}
/*
*
*/
void cJupiterModel::Run(){
    // ATJUP_FPE=1 turns the first invalid floating-point operation into a SIGFPE instead of a
    // silently propagating NaN. Run the CLI under gdb to get the exact line:
    //     cd jupiter && OMP_NUM_THREADS=1 ATJUP_FPE=1 gdb -batch -ex run -ex "bt 6"
    //         -ex "info locals" --args ../cli/jup . config_atjup.xml
    // (build with -g -O0 for line numbers). This is how the four NaN sources behind the
    // domain-wide velocity blow-up were located; NaN is invisible to printMinMax, whose
    // searchMinMax_3D compares with a bare > and so skips every non-finite cell.
    // Off by default — trapping would abort on the first harmless inf in a diagnostic field.
    if(getenv("ATJUP_FPE")) feenableexcept(FE_INVALID | FE_DIVBYZERO);

    #ifdef _OPENMP
        printf("\n\n   number of processors: %d\n\n", omp_get_num_procs());

    #pragma omp parallel
        {
            printf("   thread %d of %d in \"Desktop-Dell-XPS 8960\"\n", 
                omp_get_thread_num(), omp_get_num_threads());
        }
    #else
        printf("   OpenMP is not supported\n");
    #endif
        printf("   ended\n\n");

    int start = JupiterUtils::RunStart("JupiterCM");

    output_path = output_path + "-Jupiter";

    mkdir(output_path.c_str(), 0777);

    m_model = this;

    // ---- Metric radius (ATJUP_METRIC_RADIUS=<planet radius in km>, default off) ----
    //
    // rad.z is the radius that enters every horizontal metric factor: inv_rm, inv_rmsinthe,
    // inv_rm2sinthe2, the curvature groups of the transport terms and the horizontal parts of
    // every Laplacian. It is built as r0 = 1.0 with dr = 0.025, so it runs 1.0 to 2.0 across
    // the shell — while dr = 0.025 simultaneously means one radial step is L_atm/40 = 3.5 km.
    // The same number is therefore doing two incompatible jobs: a vertical grid spacing scaled
    // by L_atm, and a planetary radius. cJupiterModel.cpp's own declaration of r0 says as much
    // ("value much too small, Jupiter radius 72000km").
    //
    // Consequence, since Jupiter's radius over L_atm is 69911/140 = 499: every horizontal
    // derivative is about 500x larger than the geometry warrants, and every horizontal
    // Laplacian term 250000x. It bears directly on the obstacle-flank runaway — the driving
    // term w/(r sin) dw/dphi measures 20.26 there and would be 0.04 with the true radius,
    // smaller than the diffusion opposing it.
    //
    // Set the knob to the planetary radius in kilometres to give rad.z its geometric meaning,
    // r = (R + z)/L_atm. dr is untouched, so the vertical grid is exactly as before and only
    // the curvature changes. Default 0 keeps r0 = 1.0 and is bit-identical.
    //
    // NOTE this convention is inherited from ATOM, which builds rad the same way (r0 = 1.0,
    // dr = 0.025, L_atm = 400 m against an Earth radius, so the same factor ~400). Changing it
    // is a modelling decision, not a bug fix, which is why it is opt-in.
    //
    // Two places assume rad.z starts at 1 and must not be combined with this knob:
    //   - coord_stretching (default false) forms exp_rm = 1/(rm+1), which becomes 1/500 rather
    //     than 1/2 and silently rescales every radial derivative. Leave it off.
    //   - paraview_sphere_vts builds Cartesian coordinates as rad.z * (unit sphere), so the
    //     rendered shell would become a skin of relative thickness 1/500. Its only call site
    //     (FileIO_Jup.cpp) is commented out, so nothing renders wrong today; the radial, longal,
    //     zonal and panorama writers do not use rad.z at all.
    {
        static const double metric_R_km = [](){
            const char* e = getenv("ATJUP_METRIC_RADIUS"); return e ? atof(e) : 0.0; }();
        if(metric_R_km > 0.0){
            const double r0_metric = metric_R_km / L_atm;
            rad.Coordinates(im, r0_metric, dr);
            printf("      ATJUP: metric radius set from ATJUP_METRIC_RADIUS = %.0f km,"
                   " L_atm = %.1f km  ->  rad.z = %.3f .. %.3f (was 1.000 .. %.3f)\n",
                   metric_R_km, L_atm, rad.z[0], rad.z[im-1], 1.0 + (im-1)*dr);
        }
    }

    checkMetricConsistency();

    cout.precision(6);
    cout.setf(ios::fixed);

    if(!has_welcome_msg_printed)
        print_welcome_msg();

    initMsg();

    resetArrays();

    dt = 0.001;                                                         //  no dimension
    iter_n = 0;

    init_layer_heights();
    TropopauseLocation();
    init_tropopause_layers();
    VelocityInitializerJup(*this).compute();                            // construction of zonal initial velocities from measurements

    BC_Jup(*this).bcSeaMount();                                         // velocities close to surfaces, resembling a boundary layer
    computeWallViscosity();                                             // needs the SeaMount contour; fills wall_nue once
    BC_Jup(*this).bcVelSurfSur();                                       // velocities close to surfaces, resembling a boundary layer

    dampVelocities();

//    goto Printout;

    ChemistryJup(*this).ThermalPropertiesJup();

    init_temperature();

//    goto Printout;

    init_PressureStatic();

    init_vapour_cloud_ice("H2O", h2o_tropopause, 
        coeff_h2o_A, coeff_h2o_B, 
        coeff_h2o_A_i, coeff_h2o_B_i, 
        t_0_h2o, t_00_h2o, 
        ep_h2o, r_h2o, m_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o,
        h2o, h2o_cloud, h2o_ice);

//    goto Printout;

    JupiterUtils::damp_wiggles(h2o,       &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(h2o_cloud, &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(h2o_ice,   &i_topography, true, true, true);

//    goto Printout;

    init_vapour_cloud_ice("NH3", nh3_tropopause,
        coeff_nh3_A, coeff_nh3_B,
        coeff_nh3_A_i, coeff_nh3_B_i,
        t_0_nh3, t_00_nh3,
        ep_nh3, r_nh3, m_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3,
        nh3, nh3_cloud, nh3_ice);

    JupiterUtils::damp_wiggles(nh3,       &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(nh3_cloud, &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(nh3_ice,   &i_topography, true, true, true);

    init_vapour_cloud_ice("CH4", ch4_tropopause,
        coeff_ch4_A, coeff_ch4_B,
        coeff_ch4_A_i, coeff_ch4_B_i,
        t_0_ch4, t_00_ch4,
        ep_ch4, r_ch4, m_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4,
        ch4, ch4_cloud, ch4_ice);

    JupiterUtils::damp_wiggles(ch4,       &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(ch4_cloud, &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(ch4_ice,   &i_topography, true, true, true);

//    goto Printout;

    init_h2s("H2S", h2s_tropopause, 
        coeff_h2s_A, coeff_h2s_B, 
        coeff_h2s_A_i, coeff_h2s_B_i, 
        t_0_h2s, t_00_h2s, 
        ep_h2s, r_h2s, m_h2s,
        C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s, h2s);

    JupiterUtils::damp_wiggles(h2s, &i_topography, true, true, true);

//    goto Printout;

    init_nh4sh("NH4SH", nh4sh_tropopause, 
        coeff_nh4sh_A, coeff_nh4sh_B, 
        ep_nh4sh, r_nh4sh, m_nh4sh,
        C_nh4sh, L0_nh4sh, R_nh4sh, 
        del_alf_nh4sh, del_bet_nh4sh, nh4sh);

    JupiterUtils::damp_wiggles(nh4sh, &i_topography, true, true, true);

//    goto Printout;

    SaturationAdjustmentJup(*this).run("H2O",
        coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
        t_0_h2o, t_00_h2o,
        ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
        C_h2o_ice, L0_h2o_ice, del_alf_h2o_ice, del_bet_h2o_ice,
        h2o, h2o_cloud, h2o_ice);

    JupiterUtils::damp_wiggles(h2o,       &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(h2o_cloud, &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(h2o_ice,   &i_topography, true, true, true);


    SaturationAdjustmentJup(*this).run("NH3",
        coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
        t_0_nh3, t_00_nh3,
        ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
        C_nh3_ice, L0_nh3_ice, del_alf_nh3_ice, del_bet_nh3_ice,
        nh3, nh3_cloud, nh3_ice);

    JupiterUtils::damp_wiggles(nh3,       &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(nh3_cloud, &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(nh3_ice,   &i_topography, true, true, true);


    SaturationAdjustmentJup(*this).run("CH4",
        coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
        t_0_ch4, t_00_ch4,
        ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
        C_ch4_ice, L0_ch4_ice, del_alf_ch4_ice, del_bet_ch4_ice,
        ch4, ch4_cloud, ch4_ice);

    JupiterUtils::damp_wiggles(ch4,       &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(ch4_cloud, &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(ch4_ice,   &i_topography, true, true, true);


//    goto Printout;

    ChemistryJup(*this).ChemMassRateJup();
    ChemistryJup(*this).FluxLimiterNH4SH();

    JupiterUtils::damp_wiggles(massflux_h2s,   &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(massflux_nh3,   &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(massflux_nh4sh, &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(fluxlim_nh4sh,  &i_topography, true, true, true);

    ChemistryJup(*this).DiffMassFluxJup();

    JupiterUtils::damp_wiggles(difflux_h2s,     &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(difflux_nh3,     &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(difflux_nh4sh,   &i_topography, true, true, true);
    JupiterUtils::damp_wiggles(thermalmassflux, &i_topography, true, true, true);

//    goto Printout;

    Forces();
    computeMixtureDensity();
    Latent_Heat();
    init_PressureDynamic();

    // Put the horizontal wind into thermal-wind balance with the density field it has to live
    // with. Here, at the end of the initial state: the temperature, the species and p_stat are
    // final, and the boundary conditions below then tidy the edges of the field it rewrites.
    if(thermal_wind_enabled()) ThermalWindJup(*this).run();

    BC_Jup(*this).bcRadius();                                           // extrapolation in i-direction alomg grid boundaries
    BC_Jup(*this).bcTheta();                                            // extrapolation in j-direction alomg grid boundaries
    BC_Jup(*this).bcPhi();                                              // extrapolation in k-direction alomg grid boundaries

//    BC_Jup(*this).bcScalarSurfSur();                                    // scalar variable at surfaces extrapolated by von Neumann
    BC_Jup(*this).bcSolidGround();                                      // values inside mountains

    // Seed the turbulence fields from the ABL profile and prime the source terms once,
    // mirroring ATOM's TurbulenceAtm::init(). Needs the velocity field, so it comes after
    // the initial BCs.
    if(turb_enabled()) TurbulenceJup(*this).init();

    restoreVar(1.0);

//    goto Printout;

    panorama_cnt = 0;

    // One-time dump of the pristine initial state (iter=0), straight out of
    // VelocityInitializerJup, so the initial cell rotations can be inspected.
    iter_n = 0;
    writeData();

    // ---- Optional restart from a binary checkpoint ----
    // Loading here, after the full initialisation, is deliberate: the geometry (SeaMount,
    // i_topography, layer heights) and every derived constant are rebuilt from param.py as
    // usual, and only the prognostic 3D fields are then overwritten by the file. The loop
    // resumes at restart_from_iter+1, so the parity of iter_n is preserved and the
    // "even iterations do the physics" schedule continues where it left off.
    int iter_start = 1;
    if(restart_from_iter >= 0 && load_state(restart_from_iter)){
        iter_start = restart_from_iter + 1;
        restoreVar(1.0);          // refresh the n-level copies from the restored fields
    }

    for(iter_n = iter_start; iter_n <= nm; iter_n++){

        auto begin = std::chrono::high_resolution_clock::now();

        cout << endl << endl;
        cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
        cout << " 3D Jupiter iterational process" << endl;
        cout << " present state of the computation " << endl << endl 
             << " ======================== iteration number n = " << iter_n << " ========================= " << endl << endl
             << "    max total iteration number nm = " << nm << endl
             << "    checkpoint when to write 3D-panorama = " << checkpoint << endl
             << "    panorama_print = " << panorama_print << endl << endl;

        // The physics block runs on even iterations; odd iterations reuse the diagnostic
        // fields it leaves behind (massflux_*, fluxlim_nh4sh, Q_rad, Q_precip, the turbulence
        // sources, the forces). A restart deliberately stores only the PROGNOSTIC arrays, so
        // those diagnostics are still zero on the first iteration after a load — and rhs_h2s /
        // rhs_nh3 / rhs_nh4sh read massflux_* directly. Force the block on the first iteration
        // after a restart, whatever its parity, so the diagnostics are rebuilt from the
        // restored state before anything consumes them.
        if(iter_n % 2 == 0 || iter_n == iter_start){

        PressureSolverJup(*this).run();
        JupiterUtils::damp_wiggles(p_dyn, &i_topography, true, true, true);

        // Refresh the local mixture density FIRST, so the saturation adjustment, the
        // precipitation microphysics and the Lewis groups all see a rho built from the p_dyn and
        // t of this block rather than the previous one. Only matters with ATJUP_LOCAL_RHO set,
        // but the ordering should be right either way.
        computeMixtureDensity();

        if(radiation_enabled()) RadiationJup(*this).run();

        SaturationAdjustmentJup(*this).run("H2O",
            coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
            t_0_h2o, t_00_h2o,
            ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
            C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
            C_h2o_ice, L0_h2o_ice, del_alf_h2o_ice, del_bet_h2o_ice,
            h2o, h2o_cloud, h2o_ice);

        SaturationAdjustmentJup(*this).run("NH3",
            coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
            t_0_nh3, t_00_nh3,
            ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
            C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
            C_nh3_ice, L0_nh3_ice, del_alf_nh3_ice, del_bet_nh3_ice,
            nh3, nh3_cloud, nh3_ice);

        SaturationAdjustmentJup(*this).run("CH4",
            coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
            t_0_ch4, t_00_ch4,
            ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
            C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
            C_ch4_ice, L0_ch4_ice, del_alf_ch4_ice, del_bet_ch4_ice,
            ch4, ch4_cloud, ch4_ice);

        // Must stay AFTER the SaturationAdjustmentJup calls: its condensate depletion is applied
        // in place, and running it before them would let the adjustment simply undo the removal.
        if(precip_enabled()) PrecipitationJup(*this).run();   // H2O+NH3 3-cat + NH4SH settling

        // Turbulence closure. Reads the velocity field left by RK4 and the BCs, so it runs
        // after them, exactly as ATOM calls TurbulenceAtm::run() from its own iteration loop.
        if(turb_enabled()) TurbulenceJup(*this).run();

        ChemistryJup(*this).DiffMassFluxJup();                          // must precede ChemMassRateJup: massflux = w - difflux

        // Smooth difflux and thermalmassflux BEFORE difflux enters massflux assembly.
        JupiterUtils::damp_wiggles(difflux_h2s,     &i_topography, true, true, true);
        JupiterUtils::damp_wiggles(difflux_nh3,     &i_topography, true, true, true);
        JupiterUtils::damp_wiggles(difflux_nh4sh,   &i_topography, true, true, true);
        JupiterUtils::damp_wiggles(thermalmassflux, &i_topography, true, true, true);

        ChemistryJup(*this).ChemMassRateJup();
        ChemistryJup(*this).FluxLimiterNH4SH();

        // Smooth the assembled massflux.
        JupiterUtils::damp_wiggles(massflux_h2s,   &i_topography, true, true, true);
        JupiterUtils::damp_wiggles(massflux_nh3,   &i_topography, true, true, true);
        JupiterUtils::damp_wiggles(massflux_nh4sh, &i_topography, true, true, true);
        JupiterUtils::damp_wiggles(fluxlim_nh4sh,  &i_topography, true, true, true);

        Forces();
        Latent_Heat();

        }  // if loop

        RungeKuttaJup();

        BC_Jup(*this).bcRadius();                                       // extrapolation in i-direction alomg grid boundaries
        BC_Jup(*this).bcTheta();                                        // extrapolation in j-direction alomg grid boundaries
        BC_Jup(*this).bcPhi();                                          // extrapolation in k-direction alomg grid boundaries

//        BC_Jup(*this).bcScalarSurfSur();                                // scalar variable at surfaces extrapolated by von Neumann
        BC_Jup(*this).bcSolidGround();                                  // values inside mountains

        // Dry convective adjustment. It runs on the state the Runge-Kutta step just produced,
        // after the boundary conditions so that bcSolidGround has already written the solid
        // cells it must not mix across, and before restoreVar so that the adjusted temperature
        // is what the n-level copies carry into the next step. Off unless ATJUP_NONDIM is set.
        if(conv_adj_enabled()) ConvectiveAdjustmentJup(*this).run();

//        RungeKuttaJup();

        // Optional per-iteration velocity de-checkerboarding (opt-in; default off = bit-identical).
        // Use the 4th-order filter so the zonal-jet shear is preserved across iterations.
        if(shapiro_vel_inloop() > 0){
            const double s = shapiro_strength();
            JupiterUtils::damp_wiggles_ho(u, &i_topography, true, true, true, s, shapiro_vel_inloop());
            JupiterUtils::damp_wiggles_ho(v, &i_topography, true, true, true, s, shapiro_vel_inloop());
            JupiterUtils::damp_wiggles_ho(w, &i_topography, true, true, true, s, shapiro_vel_inloop());
        }

        // Floor the species at zero before the n-level copies are refreshed, so the clamped
        // values are what the next Runge-Kutta step starts from. Runs after the Shapiro pass,
        // which can itself undershoot at a sharp cloud edge. See FileIO_Jup.cpp for the
        // measurement that says a plain floor is enough here.
        clampNegativeSpecies();

        restoreVar(1.0);

        panorama_cnt++;

        if(iter_n % checkpoint == 0){
            printMinMax();
            writeData();
        }

        // Full 3D panorama .vts every 100 iterations (carries the radiation fields), independent
        // of the checkpoint/panorama_print cadence.
        if(paraview_panorama_vts_flag && iter_n % 100 == 0){
            paraview_panorama_vts(iter_n);
        }

        if(panorama_cnt == panorama_print) panorama_cnt = 1;

        // Per-iteration NaN watch (opt-in, ATJUP_NANCHECK=1): reports the first iteration at
        // which any prognostic field goes non-finite, naming the field and the cell. Costs a
        // full sweep of the restart arrays per iteration, so it is off by default.
        // Per-level momentum census (opt-in, ATJUP_WPROFILE=<stride>): every <stride>
        // iterations, print max and area-weighted mean of u,v,w for each radial level. Used to
        // decide whether the secular growth of the zonal wind is made at the radial boundary or
        // throughout the column — see momentum_profile() in FileIO_Jup.cpp.
        static const int wprof = [](){ const char* e = getenv("ATJUP_WPROFILE"); return e ? atoi(e) : 0; }();
        if(wprof > 0 && (iter_n % wprof == 0 || iter_n == iter_start)) momentum_profile(iter_n);

        static const int nan_check = [](){ const char* e = getenv("ATJUP_NANCHECK"); return e ? atoi(e) : 0; }();
        if(nan_check){
            static bool still_clean = true;
            if(still_clean && !nan_watch(iter_n)) still_clean = false;
        }

        // ---- Binary restart checkpoints ----
        // One explicit dump at checkpoint_save_iter, plus a periodic one every
        // restart_save_stride iterations. The periodic dump is written ONLY when the state is
        // clean, so a diverged run can never overwrite a good restart point — the whole value
        // of the file is that you can resume from it. Written after restoreVar so the stored
        // n-level copies are consistent with the fields they were built from.
        if(checkpoint_save_iter >= 0 && iter_n == checkpoint_save_iter)
            save_state(iter_n);

        {
            constexpr int restart_save_stride = 100;
            if(restart_save_stride > 0 && iter_n > 0 && iter_n % restart_save_stride == 0
               && iter_n != checkpoint_save_iter){
                if(restart_state_is_clean())
                    save_state(iter_n);
                else
                    cout << "      ATJUP: restart checkpoint SKIPPED at iter " << iter_n
                         << " - non-finite cell present (state not clean)" << endl;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for one time step\n", elapsed.count() * 1e-9);

    } // end iter_n

    cout << endl << "      Jupiter: run_3D_loop atm ended ..........................." << endl;



    JupiterUtils::RunEnd("JupiterCM", start);

    print_final_msg();

    return;
}
/*
*
*/
void cJupiterModel::resetArrays(){
    cout << endl << "      ATJUP: resetArrays" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    Topography.initArray_2D(jm, km, 0.0); // topography
    LatentHeat.initArray_2D(jm, km, 0.0);            // areas of higher latent heat
    Precipitation.initArray_2D(jm, km, 0.0);         // areas of higher precipitation
    precip_srf_total.initArray_2D(jm, km, 0.0);      // all-species surface precipitation [kg/m2/s]
    precip_srf_h2o.initArray_2D(jm, km, 0.0);        // H2O   rain + snow + graupel
    precip_srf_nh3.initArray_2D(jm, km, 0.0);        // NH3   rain + snow + graupel
    precip_srf_nh4sh.initArray_2D(jm, km, 0.0);      // NH4SH settling crystals
    precipitable_water.initArray_2D(jm, km, 0.0);    // areas of precipitable water in the air
    nh3_total.initArray_2D(jm, km, 0.0);             // areas of higher nh3 concentration
    nh3_cloud_total.initArray_2D(jm, km, 0.0);       // areas of higher nh3_cloud concentration
    nh3_ice_total.initArray_2D(jm, km, 0.0);         // areas of higher nh3_ice concentration
    aux_2D_v.initArray_2D(jm, km, 0.0);              // auxilliar field v
    aux_2D_w.initArray_2D(jm, km, 0.0);              // auxilliar field w
    tropopause_height.initArray_2D(jm, km, 0.0); // local height of the tropopause

    t.initArray(im, jm, km, ta);                    // temperature
    u.initArray(im, jm, km, ua);                    // u-component velocity component in r-direction
    v.initArray(im, jm, km, va);                    // v-component velocity component in theta-direction
    w.initArray(im, jm, km, wa);                    // w-component velocity component in phi-direction

    h2o.initArray(im, jm, km, 0.0);                    // water vapour
    h2o_cloud.initArray(im, jm, km, 0.0);                // cloud water
    h2o_ice.initArray(im, jm, km, 0.0);                  // cloud ice

    h2s.initArray(im, jm, km, 0.0);                    // water vapour

    nh3.initArray(im, jm, km, 0.0);                 // nh3-vapour
    nh3_cloud.initArray(im, jm, km, 0.0);           // nh3-cloud
    nh3_ice.initArray(im, jm, km, 0.0);             // nh3-ice
    ch4.initArray(im, jm, km, 0.0);                 // ch4-vapour
    ch4_cloud.initArray(im, jm, km, 0.0);           // ch4-cloud
    ch4_ice.initArray(im, jm, km, 0.0);             // ch4-ice

    nh4sh.initArray(im, jm, km, 0.0);                 // nh4sh-vapour

    tn.initArray(im, jm, km, ta);                    // temperature new
    un.initArray(im, jm, km, ua);                    // u-velocity component in r-direction new
    vn.initArray(im, jm, km, va);                    // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, wa);                    // w-velocity component in phi-direction new

    h2on.initArray(im, jm, km, 0.0);                    // water vapour new
    h2o_cloudn.initArray(im, jm, km, 0.0);                // cloud water new
    h2o_icen.initArray(im, jm, km, 0.0);                    // cloud ice new
    h2sn.initArray(im, jm, km, 0.0);                    // water vapour new
    nh3n.initArray(im, jm, km, 0.0);                // nh3 new
    nh3_cloudn.initArray(im, jm, km, 0.0);            // nh3_cloud new
    nh3_icen.initArray(im, jm, km, 0.0);            // nh3_ice new
    ch4n.initArray(im, jm, km, 0.0);                // ch4 new
    ch4_cloudn.initArray(im, jm, km, 0.0);          // ch4_cloud new
    ch4_icen.initArray(im, jm, km, 0.0);            // ch4_ice new
    nh4shn.initArray(im, jm, km, 0.0);                // nh4sh new

    fluxlim_nh4sh.initArray(im, jm, km, 0.0);  // TVD flux-limiter correction

    massflux_h2s.initArray(im, jm, km, 0.0);   // mass flux h2s
    massflux_nh3.initArray(im, jm, km, 0.0);   // mass flux nh3
    massflux_nh4sh.initArray(im, jm, km, 0.0);   // mass flux nh4sh

    difflux_h2s.initArray(im, jm, km, 0.0);   // diffusive flux h2s
    difflux_nh3.initArray(im, jm, km, 0.0);   // diffusive flux nh3
    difflux_nh4sh.initArray(im, jm, km, 0.0);   // diffusive flux nh4sh

    thermalmassflux.initArray(im, jm, km, 0.0);   // thermal massflux_h2s

    p_hydro.initArray(im, jm, km, 0.0);          // hydrostatic pressure perturbation
    p_dyn.initArray(im, jm, km, pa);                // dynamic pressure
    p_dynn.initArray(im, jm, km, pa);               // dynamic pressure (n+1)
    p_stat.initArray(im, jm, km, 1.0);                // static pressure
    rho_mix.initArray(im, jm, km, 0.0);          // local mixture density

    rho_mix.initArray(im, jm, km, 1.0);                // density of mixture

    rhs_t.initArray(im, jm, km, 0.0);                // auxilliar field RHS temperature
    rhs_u.initArray(im, jm, km, 0.0);                // auxilliar field RHS u-velocity component
    rhs_v.initArray(im, jm, km, 0.0);                // auxilliar field RHS v-velocity component
    rhs_w.initArray(im, jm, km, 0.0);                // auxilliar field RHS w-velocity component

    rhs_h2o.initArray(im, jm, km, 0.0);                // auxilliar field RHS water vapour
    rhs_h2o_cloud.initArray(im, jm, km, 0.0);            // auxilliar field RHS cloud water
    rhs_h2o_ice.initArray(im, jm, km, 0.0);                // auxilliar field RHS cloud ice
    rhs_h2s.initArray(im, jm, km, 0.0);                // auxilliar field RHS water vapour
    rhs_nh3.initArray(im, jm, km, 0.0);                // auxilliar field RHS nh3
    rhs_nh3_cloud.initArray(im, jm, km, 0.0);        // auxilliar field RHS nh3_cloud
    rhs_nh3_ice.initArray(im, jm, km, 0.0);            // auxilliar field RHS nh3_ice
    rhs_ch4.initArray(im, jm, km, 0.0);                // auxilliar field RHS ch4
    rhs_ch4_cloud.initArray(im, jm, km, 0.0);        // auxilliar field RHS ch4_cloud
    rhs_ch4_ice.initArray(im, jm, km, 0.0);            // auxilliar field RHS ch4_ice
    rhs_nh4sh.initArray(im, jm, km, 0.0);                // auxilliar field RHS nh4sh
    rhs_tke.initArray(im, jm, km, 0.0);                // auxilliar field RHS turbulent kinetic energy
    rhs_dis.initArray(im, jm, km, 0.0);                // auxilliar field RHS dissipation

    aux.initArray(im, jm, km, 0.0);                // auxilliar field u-velocity component
    aux_u.initArray(im, jm, km, 0.0);                // auxilliar field u-velocity component
    aux_v.initArray(im, jm, km, 0.0);                // auxilliar field v-velocity component
    aux_w.initArray(im, jm, km, 0.0);                // auxilliar field w-velocity component

    Q_Latent.initArray(im, jm, km, 0.0);                // latent heat
    Q_Sensible.initArray(im, jm, km, 0.0);            // sensible heat
    Q_rad.initArray(im, jm, km, 0.0);                // radiative heating rate [W/m3]
    radiation.initArray(im, jm, km, 0.0);            // layer-centre net radiative flux [W/m2]
    epsilon.initArray(im, jm, km, 0.0);                // layer emissivity
    P_rain.initArray(im, jm, km, 0.0);               // H2O rain    precipitation flux [kg/m2/s]
    P_snow.initArray(im, jm, km, 0.0);               // H2O snow    precipitation flux [kg/m2/s]
    P_graupel.initArray(im, jm, km, 0.0);            // H2O graupel precipitation flux [kg/m2/s]
    P_nh3_rain.initArray(im, jm, km, 0.0);           // NH3 rain    precipitation flux [kg/m2/s]
    P_nh3_snow.initArray(im, jm, km, 0.0);           // NH3 snow    precipitation flux [kg/m2/s]
    P_nh3_graupel.initArray(im, jm, km, 0.0);        // NH3 graupel precipitation flux [kg/m2/s]
    P_nh4sh.initArray(im, jm, km, 0.0);              // NH4SH crystal sedimentation flux [kg/m2/s]
    Q_precip.initArray(im, jm, km, 0.0);             // latent heating rate from precip [W/m3]
    tke.initArray(im, jm, km, 0.0);                   // turbulent kinetic energy k*
    tken.initArray(im, jm, km, 0.0);                  // k* at time level n
    dis.initArray(im, jm, km, 1.0e-10);               // epsilon* or omega*
    disn.initArray(im, jm, km, 1.0e-10);              // dis at time level n
    nue.initArray(im, jm, km, 0.0);                   // eddy viscosity nue*
    prod.initArray(im, jm, km, 0.0);                  // production contraction P_k
    tke_source.initArray(im, jm, km, 0.0);            // k source term
    dis_source.initArray(im, jm, km, 0.0);            // dis source term
    vel_star.initArray_2D(jm, km, 0.0);               // friction velocity u_tau [m/s]
    CoriolisForce.initArray(im, jm, km, 0.0);        // Coriolis force
    CentrifugalForce.initArray(im, jm, km, 0.0);             // centrifugal force
    BuoyancyForce.initArray(im, jm, km, 0.0);        // buoyancy force, Boussinesque approximation
    PresGradForce.initArray(im, jm, km, 0.0);// pressure gradient force
    SeaMount.initArray(im, jm, km, 0.0);             // sea mount contour

    w_nh3.initArray(im, jm, km, 0.0);                // chemical reaction rate of nh3
    w_h2s.initArray(im, jm, km, 0.0);                // chemical reaction rate of h2s
    w_nh4sh.initArray(im, jm, km, 0.0);                // chemical reaction rate of nh4sh

    j_nh3.initArray(im, jm, km, 0.0);                // ordinary-diffusion mass flux of nh3
    j_h2s.initArray(im, jm, km, 0.0);                // ordinary-diffusion mass flux of h2s
    j_nh4sh.initArray(im, jm, km, 0.0);              // ordinary-diffusion mass flux of nh4sh

    jT_nh3.initArray(im, jm, km, 0.0);               // thermo-diffusion mass flux of nh3
    jT_h2s.initArray(im, jm, km, 0.0);               // thermo-diffusion mass flux of h2s
    jT_nh4sh.initArray(im, jm, km, 0.0);             // thermo-diffusion mass flux of nh4sh

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for resetArrays\n", elapsed.count() * 1e-9);

    cout << "      ATJUP: resetArrays ended" << endl;
    return;
}
/*
*
*/
void cJupiterModel::restoreVar(double coeff){
//    cout << endl << "      ATJUP: restoreVar" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                tn.x[i][j][k] = coeff * t.x[i][j][k];
                un.x[i][j][k] = coeff * u.x[i][j][k];
                vn.x[i][j][k] = coeff * v.x[i][j][k];
                wn.x[i][j][k] = coeff * w.x[i][j][k];
                h2on.x[i][j][k] = coeff * h2o.x[i][j][k];
                h2o_cloudn.x[i][j][k] = coeff * h2o_cloud.x[i][j][k];
                h2o_icen.x[i][j][k] = coeff * h2o_ice.x[i][j][k];
                h2sn.x[i][j][k] = coeff * h2s.x[i][j][k];
                nh3n.x[i][j][k] = coeff * nh3.x[i][j][k];
                nh3_cloudn.x[i][j][k] = coeff * nh3_cloud.x[i][j][k];
                nh3_icen.x[i][j][k] = coeff * nh3_ice.x[i][j][k];
                ch4n.x[i][j][k] = coeff * ch4.x[i][j][k];
                ch4_cloudn.x[i][j][k] = coeff * ch4_cloud.x[i][j][k];
                ch4_icen.x[i][j][k] = coeff * ch4_ice.x[i][j][k];
                nh4shn.x[i][j][k] = coeff * nh4sh.x[i][j][k];
                // k* and dis* are now prognostic in RungeKuttaJup (rhs_tke / rhs_dis), so their
                // time-level-n copies have to be refreshed here exactly like every other
                // transported variable — this is ATOM's UtilsAtm::storeIntermediateData3D.
                tken.x[i][j][k] = coeff * tke.x[i][j][k];
                disn.x[i][j][k] = coeff * dis.x[i][j][k];
            }
        }
    }

//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
//    printf(" time measured: %.3f seconds for restoreVar\n", elapsed.count() * 1e-9);

//    cout << "      ATJUP: restoreVar ended" << endl;

    return;
}
/*
*
*/
/*
    cout << endl << "      ATJUP: thermodynamic fft_gaussian_filter_3d in Run begin ......................." << endl;

    fft_gaussian_filter_3d(nh3,1);
    fft_gaussian_filter_3d(nh3_cloud,1);
    fft_gaussian_filter_3d(nh3_ice,1);

    fft_gaussian_filter_3d(nh4sh,1);
//    fft_gaussian_filter_3d(nh4sh_cloud,1);

    fft_gaussian_filter_3d(h2s,1);
    fft_gaussian_filter_3d(h2s_cloud,1);
    fft_gaussian_filter_3d(h2s_ice,1);

    fft_gaussian_filter_3d(h2o,1);
    fft_gaussian_filter_3d(h2o_cloud,1);
    fft_gaussian_filter_3d(h2o_ice,1);

    cout << endl << "      ATJUP: thermodynamic fft_gaussian_filter_3d in Run end ......................." << endl;
*/


