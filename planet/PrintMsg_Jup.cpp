#include "cJupiterModel.h"
#include "Reporting.h"

using namespace std;
using namespace JupiterUtils;

/*
* The min/max report machinery and steadyQuery now live in the SHARED Reporting.h, as
* Reporting<Planet>. What stayed here is what is genuinely Jupiter's: which fields printMinMax
* lists, in what unit, and the welcome/final text. See Reporting.h for why the split falls there
* — in particular why printMinMax itself is NOT shared, the two models having settled the species
* units differently.
*/
// UNITS OF THE SPECIES FIELDS — settled here, because the model held three readings at once.
//
// c is a mass DENSITY in kg/m3. init_vapour_cloud_ice builds it as magnus*r_mix*q_sat and caps
// it at r_max, and the r_* constants are documented as "density ... in kg/m3"; PrecipitationJup
// states it outright ("q is already a mass concentration [kg/m3]"); the NH4SH Stokes flux in
// rhs_nh4sh assumes the same. Against that, this block multiplied every species by r_mix again
// before printing "kg/m3", which made the displayed value 1.2844x too large — 28 %. The
// coefficient is now 1.0 (and 1e3 where grams are wanted), so the printout finally shows the
// field the model actually carries. RadiationJup was the third reading and is corrected there.
void cJupiterModel::printMinMax(){
    // The 39 fields every model shares are the SHARED Reporting<Planet>; see
    // print_minmax_common() there for why they could not be shared until the units were
    // settled. Below is only what is ATJUP's own.
    Reporting<cJupiterModel>(*this).print_minmax_common();

    // t*t_ref is the physical temperature in K (t_ref = 165, the scale SaturationAdjustmentJup
    // and RadiationJup compare against), so degC needs -273.15. This used to subtract 165.0,
    // i.e. it printed the offset from t_ref while labelling it degC — every temperature in the
    // log read ~108 K too high. Solid/masked cells (t = 0) now show as -273.15 degC.
    cout << endl;

    cout << endl;

    // p_dyn is stored as the nondimensional kinematic pressure; p_dyn_to_bar() makes the
    // printed number match its unit label. See the note in cJupiterModel.h.
    cout << endl;

    cout << endl;

    // CH4 was added to the transport, chemistry and saturation pipeline but never to this
    // report, so the only place its fields showed up was ParaView. NOTE for reading the
    // numbers: ch4_cloud and ch4_ice are expected to stay ZERO on Jupiter. Methane condenses
    // near 80-90 K at these pressures and Jupiter's coldest level is ~110 K, so CH4 is a
    // well-mixed non-condensing gas here (it condenses on Uranus and Neptune, not Jupiter) —
    // hence the "NO saturation found in SaturationAdjustment of CH4" line each iteration.
    cout << endl;

    cout << endl;

    cout << endl;

    cout << endl;

    cout << endl;

    cout << endl << " Precipitation " << endl;
    searchMinMax_3D(" max 3D P_rain ", " min 3D P_rain ", " kg/m2/s", P_rain, 1.0);
    searchMinMax_3D(" max 3D P_snow ", " min 3D P_snow ", " kg/m2/s", P_snow, 1.0);
    searchMinMax_3D(" max 3D P_graupel ", " min 3D P_graupel ", " kg/m2/s", P_graupel, 1.0);
    searchMinMax_3D(" max 3D P_nh3_rain ", " min 3D P_nh3_rain ", " kg/m2/s", P_nh3_rain, 1.0);
    searchMinMax_3D(" max 3D P_nh3_snow ", " min 3D P_nh3_snow ", " kg/m2/s", P_nh3_snow, 1.0);
    searchMinMax_3D(" max 3D P_nh4sh ", " min 3D P_nh4sh ", " mg/m3", P_nh4sh, 1e6);
    searchMinMax_3D(" max 3D Q_precip ", " min 3D Q_precip ", " W/m3", Q_precip, 1.0);
    // All-species surface map, in mm/day so it can be compared against the ~1 mm/d that
    // Jupiter's energy budget allows (see the calibration note in PrecipitationJup.h).
    searchMinMax_2D(" max 2D precip surface total ", " min 2D precip surface total ", " mm/d", precip_srf_total, 86400.0);
    searchMinMax_2D(" max 2D precip surface H2O ",   " min 2D precip surface H2O ",   " mm/d", precip_srf_h2o,   86400.0);
    searchMinMax_2D(" max 2D precip surface NH3 ",   " min 2D precip surface NH3 ",   " mm/d", precip_srf_nh3,   86400.0);
    searchMinMax_2D(" max 2D precip surface NH4SH ", " min 2D precip surface NH4SH ", " mm/d", precip_srf_nh4sh, 86400.0);
    cout << endl;

    cout << endl << " Radiation " << endl;
    searchMinMax_3D(" max 3D net radiation ", " min 3D net radiation ", " W/m2", radiation, 1.0);
    searchMinMax_3D(" max 3D Q_rad ", " min 3D Q_rad ", " W/m3", Q_rad, 1.0);
    searchMinMax_3D(" max 3D emissivity ", " min 3D emissivity ", " /", epsilon, 1.0);

    // k-omega SST turbulence (physical units; zero unless ATJUP_TURB is set)
    cout << endl << " Turbulence (dimensionless: k*, eps*/omega*, nue*) " << endl;
    searchMinMax_3D(" max 3D tke ", " min 3D tke ", " /", tke, 1.0);
    searchMinMax_3D(" max 3D dis ", " min 3D dis ", " /", dis, 1.0);
    searchMinMax_3D(" max 3D nue ", " min 3D nue ", " /", nue, 1.0);
    searchMinMax_3D(" max 3D prod ", " min 3D prod ", " /", prod, 1.0);
    searchMinMax_3D(" max 3D tke_source ", " min 3D tke_source ", " /", tke_source, 1.0);
    searchMinMax_3D(" max 3D dis_source ", " min 3D dis_source ", " /", dis_source, 1.0);
    searchMinMax_2D(" max 2D vel_star ", " min 2D vel_star ", " m/s", vel_star, 1.0);

    // Equatorial column profile (j=jm/2, k=km/2), top -> bottom, for a direct
    // check of the radiation / Q_rad fields against the actual T(p).
    {
        const int j0 = jm / 2, k0 = km / 2;
        cout << endl << " Equatorial column  (j=" << j0 << ", k=" << k0
             << ")   top -> bottom" << endl;
        printf("   %3s  %10s  %8s  %8s  %12s  %14s\n",
               "i", "p[bar]", "T[K]", "eps", "netRad[W/m2]", "Q_rad[W/m3]");
        for(int i = im - 1; i >= 0; i--){
            printf("   %3d  %10.4f  %8.2f  %8.4f  %12.4f  %14.4e\n",
                   i, p_stat.x[i][j0][k0], t.x[i][j0][k0] * t_ref,
                   epsilon.x[i][j0][k0], radiation.x[i][j0][k0], Q_rad.x[i][j0][k0]);
        }
    }
    cout << endl << endl;

    reportClampBudget();
}
/*
*
*/
/*
*
*/
// Forwarders to the shared Reporting<cJupiterModel>. The bodies used to be here in full and were
// the highest-overlap pair between the two models (83.8 % line-for-line).
void cJupiterModel::searchMinMax_3D(string name_maxValue, string name_minValue,
    string name_unitValue, Array &value_D, double coeff,
    std::function< double(double) > lambda, bool print_heading){
    Reporting<cJupiterModel>(*this).searchMinMax_3D(name_maxValue, name_minValue,
        name_unitValue, value_D, coeff, lambda, print_heading);
}
/*
*
*/
void cJupiterModel::searchMinMax_2D(string name_maxValue, string name_minValue,
    string name_unitValue, Array_2D &value, double coeff){
    Reporting<cJupiterModel>(*this).searchMinMax_2D(name_maxValue, name_minValue,
        name_unitValue, value, coeff);
}
void cJupiterModel::print_welcome_msg(){
    if(verbose){
        cout << endl << endl << endl;
        cout << "***** Atmosphere Jupiter General Circulation Model(ATJUP) applied to laminar flow" << endl;
        cout << "***** program for the computation of Jupiter-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 6 additional transport equations to describe the water vapour, cloud water, cloud ice and nh3 vapour, nh3 cloud and nh3 ice" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonally constant" << endl;
        cout << "***** water and nh3 vapour distribution given by Clausius-Claperon equation for the partial pressure" << endl;
        cout << "***** water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model" << endl;
        cout << "***** two category ice scheme for cold clouds applying parameterization schemes provided by the COSMO code(German Weather Forecast)" << endl;
        cout << "***** rain and snow precipitation solved by column equilibrium applying the diagnostic equations" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)" << endl << endl;
        cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
        has_welcome_msg_printed = true;
    }
    return;
}
/*
*
*/
void cJupiterModel::initMsg(){
    cout << "  present state of the computation " << endl << "  current number of iterations iter_n = " << iter_n << endl << endl;
    return;
}
/*
*
*/
void cJupiterModel::print_final_msg(){
    cout << endl << "***** end of the JupiterAtmosphere General Circulation Modell(ATJUP) *****" << endl << endl;
    if(n == nm)   cout <<  "***** number of maximum number of time steps reached     nm = " << nm << endl << endl;
}
/*
*
*/
// See Reporting.h for what this reports and the three defects that were fixed when it was
// revived. It must run BEFORE restoreVar(), or every number it prints is identically zero.
void cJupiterModel::steadyQuery(){
    Reporting<cJupiterModel>(*this).steadyQuery();
}
