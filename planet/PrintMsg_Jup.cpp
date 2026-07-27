#include "cJupiterModel.h"

using namespace std;
using namespace JupiterUtils;

namespace{
    string heading_1 = " printout of maximum and minimum values of properties at their locations: latitude, longitude, level";
    string heading_2 = " results based on three dimensional considerations of the problem";
    string level = "km";
}
/*
*
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

    cout << endl << endl << " Courant time step   dt = " << dt << endl << endl;



    cout << endl << endl << " Temperatures " << endl;
    // t*t_ref is the physical temperature in K (t_ref = 165, the scale SaturationAdjustmentJup
    // and RadiationJup compare against), so degC needs -273.15. This used to subtract 165.0,
    // i.e. it printed the offset from t_ref while labelling it degC — every temperature in the
    // log read ~108 K too high. Solid/masked cells (t = 0) now show as -273.15 degC.
    searchMinMax_3D(" max 3D temperature ", " min 3D temperature ","  degC",
        t, t_ref, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max 3D thermalflux ", " min 3D thermalflux ", " kJ/(m5K)", thermalmassflux, 1e-3);
    cout << endl;

    cout << endl << " Velocities " << endl;
    searchMinMax_3D(" max 3D u-component ", " min 3D u-component ", "m/s", u, u_0);
    searchMinMax_3D(" max 3D v-component ", " min 3D v-component ", "m/s", v, u_0);
    searchMinMax_3D(" max 3D w-component ", " min 3D w-component ", "m/s", w, u_0);
    cout << endl;

    cout << endl << " Pressures " << endl;
    searchMinMax_3D(" max 3D pressure dynamic ", " min 3D pressure dynamic ", "bar", p_dyn, 1.0);
    searchMinMax_3D(" max 3D pressure static ", " min 3D pressure static ", "bar", p_stat, 1.0);
    cout << endl;

    cout << endl << " Water " << endl;
    searchMinMax_3D(" max 3D h2o ",  " min 3D h2o ", " kg/m3", h2o, 1.0);
    searchMinMax_3D(" max 3D h2o_cloud ", " min 3D h2o_cloud ", " kg/m3", h2o_cloud, 1.0);
    searchMinMax_3D(" max 3D h2o_ice ", " min 3D h2o_ice ", " kg/m3", h2o_ice, 1.0);
    cout << endl;

    // CH4 was added to the transport, chemistry and saturation pipeline but never to this
    // report, so the only place its fields showed up was ParaView. NOTE for reading the
    // numbers: ch4_cloud and ch4_ice are expected to stay ZERO on Jupiter. Methane condenses
    // near 80-90 K at these pressures and Jupiter's coldest level is ~110 K, so CH4 is a
    // well-mixed non-condensing gas here (it condenses on Uranus and Neptune, not Jupiter) —
    // hence the "NO saturation found in SaturationAdjustment of CH4" line each iteration.
    cout << endl << " Methane " << endl;
    searchMinMax_3D(" max 3D ch4 ", " min 3D ch4 ", " kg/m3", ch4, 1.0);
    searchMinMax_3D(" max 3D ch4_cloud ", " min 3D ch4_cloud ", " kg/m3", ch4_cloud, 1.0);
    searchMinMax_3D(" max 3D ch4_ice ", " min 3D ch4_ice ", " kg/m3", ch4_ice, 1.0);
    cout << endl;

    cout << endl << " Hydrogen Sulfide " << endl;
    searchMinMax_3D(" max 3D h2s ",  " min 3D h2s ", " kg/m3", h2s, 1.0);
    searchMinMax_3D(" max 3D w_h2s ", " min 3D w_h2s ", " kg/m3s", w_h2s, 1.0);
    searchMinMax_3D(" max 3D j_h2s ", " min 3D j_h2s ", " kg/m4", j_h2s, 1.0);
    searchMinMax_3D(" max 3D jT_h2s ", " min 3D jT_h2s ", " kg/m4", jT_h2s, 1.0);
    searchMinMax_3D(" max 3D massflux_h2s ", " min 3D massflux_h2s ", " kg/m3s", massflux_h2s, 1.0);
    searchMinMax_3D(" max 3D diff_h2s ", " min 3D diff_h2s ", " kg/m3s", difflux_h2s, 1.0);
    cout << endl;

    cout << endl << " Ammonia " << endl;
    searchMinMax_3D(" max 3D nh3 ",  " min 3D nh3 ", " kg/m3", nh3, 1.0);
    searchMinMax_3D(" max 3D nh3_cloud ", " min 3D nh3_cloud ", " kg/m3", nh3_cloud, 1.0);
    searchMinMax_3D(" max 3D nh3_ice ", " min 3D nh3_ice ", " kg/m3", nh3_ice, 1.0);
    searchMinMax_3D(" max 3D w_nh3 ", " min 3D w_nh3 ", " kg/m3s", w_nh3, 1.0);
    searchMinMax_3D(" max 3D j_nh3 ", " min 3D j_nh3 ", " kg/m4", j_nh3, 1.0);
    searchMinMax_3D(" max 3D jT_nh3 ", " min 3D jT_nh3 ", " kg/m4", jT_nh3, 1.0);
    searchMinMax_3D(" max 3D massflux_nh3 ", " min 3D massflux_nh3 ", " kg/m3s", massflux_nh3, 1.0);
    searchMinMax_3D(" max 3D diff_nh3 ", " min 3D diff_nh3 ", " kg/(m3s)", difflux_nh3, 1.0);
    cout << endl;

    cout << endl << " Ammonia Hydrosulfide " << endl;
    searchMinMax_3D(" max 3D nh4sh ",  " min 3D nh4sh ", " g/m3", nh4sh, 1e3);
    searchMinMax_3D(" max 3D w_nh4sh ", " min 3D w_nh4sh ", " g/m3s", w_nh4sh, 1e3);
    searchMinMax_3D(" max 3D massflux_nh4sh ", " min 3D massflux_nh4sh ", " g/m3s", massflux_nh4sh, 1e3);
    searchMinMax_3D(" max 3D j_nh4sh ", " min 3D j_nh4sh ", " g/m4", j_nh4sh, 1e3);
    searchMinMax_3D(" max 3D jT_nh4sh ", " min 3D jT_nh4sh ", " g/m4", jT_nh4sh, 1e3);
    searchMinMax_3D(" max 3D diff_nh4sh ", " min 3D diff_nh4sh ", " g/(m3s)", difflux_nh4sh, 1e3);
    cout << endl;

    cout << endl << " Forces " << endl;
    searchMinMax_3D(" max 3D Coriolis force ", " min 3D Coriolis force ", " mN/m3", CoriolisForce, 1e3);
    searchMinMax_3D(" max 3D centrifugal force ", " min 3D centrifugal force ", " mN/m3", CentrifugalForce, 1e3);
    searchMinMax_3D(" max 3D buoyancy force ", " min 3D buoyancy force ", " N/m3", BuoyancyForce, 1.0);
    searchMinMax_3D(" max 3D presgrad force ", " min 3D presgrad force ", " N/m3", PresGradForce, 1.0);
    cout << endl;

    cout << endl << " Energies " << endl;
    searchMinMax_3D(" max 3D sensible heat ", " min 3D sensible heat ", " W/m3", Q_Sensible, 1.0);
    searchMinMax_3D(" max 3D latent heat ", " min 3D latent heat ", " W/m3", Q_Latent, 1.0);
    cout << endl;

    cout << endl << " Precipitation " << endl;
    searchMinMax_3D(" max 3D P_rain ", " min 3D P_rain ", " kg/m2/s", P_rain, 1.0);
    searchMinMax_3D(" max 3D P_snow ", " min 3D P_snow ", " kg/m2/s", P_snow, 1.0);
    searchMinMax_3D(" max 3D P_graupel ", " min 3D P_graupel ", " kg/m2/s", P_graupel, 1.0);
    searchMinMax_3D(" max 3D P_nh3_rain ", " min 3D P_nh3_rain ", " kg/m2/s", P_nh3_rain, 1.0);
    searchMinMax_3D(" max 3D P_nh3_snow ", " min 3D P_nh3_snow ", " kg/m2/s", P_nh3_snow, 1.0);
    searchMinMax_3D(" max 3D P_nh4sh ", " min 3D P_nh4sh ", " kg/m2/s", P_nh4sh, 1.0);
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
void cJupiterModel::searchMinMax_3D(string name_maxValue, string name_minValue, 
    string name_unitValue, Array &value_D, double coeff, 
    std::function< double(double) > lambda, bool print_heading){
    double maxValue = value_D.x[0][0][0];
    double minValue = value_D.x[0][0][0];
    int imax = 0;
    int jmax = 0;
    int kmax = 0;
    int imin = 0;
    int jmin = 0;
    int kmin = 0;  
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(value_D.x[i][j][k] > maxValue){
                    maxValue = value_D.x[i][j][k];
                    imax = i;
                    jmax = j;
                    kmax = k;
                }else if(value_D.x[i][j][k] < minValue){
                    minValue = value_D.x[i][j][k];
                    imin = i;
                    jmin = j;
                    kmin = k;
                }
            }
        }
    }
    int imax_level = imax * (int)L_atm/(im-1);
    int imin_level = imin * (int)L_atm/(im-1);
    //  maximum latitude and longitude units recalculated
    HemisphereCoords coords = convert_coords(kmax, jmax);
    int jmax_deg = coords.lat;
    string deg_lat_max = coords.north_or_south;
    int kmax_deg = coords.lon;
    string deg_lon_max = coords.east_or_west;
    //  minimum latitude and longitude units recalculated
    coords = convert_coords(kmin, jmin);
    int jmin_deg = coords.lat;
    string deg_lat_min= coords.north_or_south;
    int kmin_deg = coords.lon;
    string deg_lon_min = coords.east_or_west;
    cout.precision(6);
    if(print_heading){
        cout << endl << heading_1 << endl << heading_2 << endl << endl;
    }
    maxValue = lambda(maxValue * coeff);
    minValue = lambda(minValue * coeff);
    cout << setiosflags(ios::left) << setw(26) << setfill('.') << name_maxValue << " = " <<
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << maxValue << setw(12) <<
        name_unitValue << setw(5) << jmax_deg << setw(3) << deg_lat_max << setw(4) << kmax_deg <<
        setw(3) << deg_lon_max << setw(6) << imax_level << setw(2) << level << "   " <<
        setiosflags(ios::left) << setw(26) << setfill('.') << name_minValue << " = "<<
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << minValue << setw(12) <<
        name_unitValue << setw(5)  << jmin_deg << setw(3) << deg_lat_min << setw(4) << kmin_deg <<
        setw(3) << deg_lon_min  << setw(6) << imin_level << setw(2) << level << endl;
}
/*
*
*/
void cJupiterModel::searchMinMax_2D(string name_maxValue, string name_minValue, 
    string name_unitValue, Array_2D &value, double coeff){
    double minValue = value.y[0][0];
    double maxValue = value.y[0][0];
    int jmax = 0;
    int kmax = 0;
    int jmin = 0;
    int kmin = 0;  
    for(int j = 1; j < jm-1; j++){
        for(int k = 1; k < km-1; k++){
            if(value.y[j][k] > maxValue){
                maxValue = value.y[j][k];
                jmax = j;
                kmax = k;
            }else if(value.y[j][k] < minValue){
                minValue = value.y[j][k];
                jmin = j;
                kmin = k;
            }
        }
    }
    int imax_level = 0;
    int imin_level = 0;
    //  maximum latitude and longitude units recalculated
    HemisphereCoords coords = convert_coords(kmax, jmax);
    int jmax_deg = coords.lat;
    string deg_lat_max = coords.north_or_south;
    int kmax_deg = coords.lon;
    string deg_lon_max = coords.east_or_west;
    //  minimum latitude and longitude units recalculated
    coords = convert_coords(kmin, jmin);
    int jmin_deg = coords.lat;
    string deg_lat_min= coords.north_or_south;
    int kmin_deg = coords.lon;
    string deg_lon_min = coords.east_or_west;
    cout.precision(6);
    maxValue = maxValue * coeff;
    minValue = minValue * coeff;
    cout << setiosflags(ios::left) << setw(26) << setfill('.') << name_maxValue << " = " <<
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << maxValue << setw(12) <<
        name_unitValue << setw(5) << jmax_deg << setw(3) << deg_lat_max << setw(4) << kmax_deg <<
        setw(3) << deg_lon_max << setw(6) << imax_level << setw(2) << level << "   " <<
        setiosflags(ios::left) << setw(26) << setfill('.') << name_minValue << " = "<<
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << minValue << setw(12) <<
        name_unitValue << setw(5)  << jmin_deg << setw(3) << deg_lat_min << setw(4) << kmin_deg <<
        setw(3) << deg_lon_min  << setw(6) << imin_level << setw(2) << level << endl;
}
/*
*
*/
double cJupiterModel::out_maxValue() const{
    return maxValue;
}
/*
*
*/
double cJupiterModel::out_minValue() const{
    return minValue;
}
/*
*
*/
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
void cJupiterModel::steadyQuery(){
    int i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, 
        j_c, k_c, i_cloud, j_cloud, k_cloud, i_ice, j_ice, k_ice, i_nh3, 
        j_nh3, k_nh3, i_nh3_cloud, j_nh3_cloud, k_nh3_cloud, i_nh3_ice, 
        j_nh3_ice, k_nh3_ice, i_nh4sh, j_nh4sh, k_nh4sh, i_p, j_p, k_p;
    int i_loc, j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg;
    double max_u, max_v, max_w, max_t, max_c, max_cloud, max_ice, max_p, 
        max_nh3, max_nh3_cloud, max_nh3_ice, max_nh4sh;
    double min_u, min_v, min_w, min_t, min_c, min_cloud, min_ice, 
        min_nh3, min_nh3_cloud, min_nh3_ice, min_nh4sh, min_p;
    double Value;
    string name_Value;
    string level, deg_north, deg_south, deg_west, deg_east, deg_lat, 
        deg_lon, heading;
    min_u = max_u = 0.;
    min_v = max_v = 0.;
    min_w = max_w = 0.;
    min_t = max_t = 0.;
    min_c = max_c = 0.;
    min_p = max_p = 0.;
    min_cloud = max_cloud = 0.;
    min_ice = max_ice = 0.;
    min_nh3 = max_nh3 = 0.;
    min_nh3_cloud = max_nh3_cloud = 0.;
    min_nh3_ice = max_nh3_ice = 0.;
    double sinthe = 0., costhe = 0., rmsinthe = 0.;
    double dudr = 0., dvdthe = 0., dwdphi = 0.;
    double residuum = 0.;
    double minimum = 0.;
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            if(costhe_abs() && j > 90) costhe = - costhe;
            rmsinthe = rad.z[i] * sinthe;
            for(int k = 1; k < km-1; k++){
                dudr = (u.x[i+1][j][k] - u.x[i-1][j][k])/(2. * dr);
                dvdthe = (v.x[i][j+1][k] - v.x[i][j-1][k])/(2. * dthe);
                dwdphi = (w.x[i][j][k+1] - w.x[i][j][k-1])/(2. * dphi);
                residuum = dudr + 2. * u.x[i][j][k]/rad.z[i] + dvdthe/rad.z[i]
                    + costhe/rmsinthe * v.x[i][j][k] + dwdphi/rmsinthe;
                if(fabs(residuum) >= minimum){
                    minimum = residuum;
                    i_res = i;
                    j_res = j;
                    k_res = k;
                }
            }
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                p_dynn.x[i][j][k] = p_dyn.x[i][j][k];
                max_p = fabs(p_dyn.x[i][j][k] - p_dynn.x[i][j][k]);
                if(max_p >= min_p){
                    min_p = max_p;
                    i_p = i;
                    j_p = j;
                    k_p = k;
                }
                max_u = fabs(u.x[i][j][k] - un.x[i][j][k]);
                if(max_u >= min_u){
                    min_u = max_u;
                    i_u = i;
                    j_u = j;
                    k_u = k;
                }
                max_v = fabs(v.x[i][j][k] - vn.x[i][j][k]);
                if(max_v >= min_v){
                    min_v = max_v;
                    i_v = i;
                    j_v = j;
                    k_v = k;
                }
                max_w = fabs(w.x[i][j][k] - wn.x[i][j][k]);
                if(max_w >= min_w){
                    min_w = max_w;
                    i_w = i;
                    j_w = j;
                    k_w = k;
                }
                max_t = fabs(t.x[i][j][k] - tn.x[i][j][k]);
                if(max_t >= min_t){
                    min_t = max_t;
                    i_t = i;
                    j_t = j;
                    k_t = k;
                }
                max_c = fabs(h2o.x[i][j][k] - h2on.x[i][j][k]);
                if(max_c >= min_c){
                    min_c = max_c;
                    i_c = i;
                    j_c = j;
                    k_c = k;
                }
                max_cloud = fabs(h2o_cloud.x[i][j][k] - h2o_cloudn.x[i][j][k]);
                if(max_cloud >= min_cloud){
                    min_cloud = max_cloud;
                    i_cloud = i;
                    j_cloud = j;
                    k_cloud = k;
                }
                max_ice = fabs(h2o_ice.x[i][j][k] - h2o_icen.x[i][j][k]);
                if(max_ice >= min_ice){
                    min_ice = max_ice;
                    i_ice = i;
                    j_ice = j;
                    k_ice = k;
                }
                max_nh3 = fabs(nh3.x[i][j][k] - nh3n.x[i][j][k]);
                if(max_nh3 >= min_nh3){
                    min_nh3= max_nh3;
                    i_nh3 = i;
                    j_nh3 = j;
                    k_nh3 = k;
                }
                max_nh3_cloud = fabs(nh3_cloud.x[i][j][k] - nh3_cloudn.x[i][j][k]);
                if(max_nh3_cloud >= min_nh3_cloud){
                    min_nh3_cloud= max_nh3_cloud;
                    i_nh3_cloud = i;
                    j_nh3_cloud = j;
                    k_nh3_cloud = k;
                }
                max_nh3_ice = fabs(nh3_ice.x[i][j][k] - nh3_icen.x[i][j][k]);
                if(max_nh3_ice >= min_nh3_ice){
                    min_nh3_ice= max_nh3_ice;
                    i_nh3_ice = i;
                    j_nh3_ice = j;
                    k_nh3_ice = k;
                }
                max_nh4sh = fabs(nh4sh.x[i][j][k] - nh4shn.x[i][j][k]);
                if(max_nh4sh >= min_nh4sh){
                    min_nh4sh= max_nh4sh;
                    i_nh4sh = i;
                    j_nh4sh = j;
                    k_nh4sh = k;
                }
            }
        }
    }
cout << "   " << i_p << "   " << j_p << "   " << k_p << "   " << max_p << "   " << min_p << endl;
cout << "   " << i_t << "   " << j_t << "   " << k_t << "   " << max_t << "   " << min_t << endl;
    cout.precision(6);
    cout.setf(ios::fixed);
    cout << endl << endl;
    cout << "      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    cout << "      3D ATJUP iterational process" << endl;
    cout << "      max total iteration number nm = " << nm << endl;
    cout << "      n = " << n << endl << endl;
    cout << endl;
    heading = " 3D iterational process ";
    cout << endl << endl << heading << endl << endl;
    level = "km";
    deg_north = "°N";
    deg_south = "°S";
    deg_west = "°W";
    deg_east = "°E";
    int choice = {1};
    preparation:
    switch(choice){
        case 1 :    name_Value = " residuum: continuity equation ";
                        Value = minimum;
                        i_loc = i_res;
                        j_loc = j_res;
                        k_loc = k_res;
                        break;
        case 2 :    name_Value = " dp: pressure Poisson equation ";

                        Value = min_p;
                        i_loc = i_p;
                        j_loc = j_p;
                        k_loc = k_p;
                        break;
        case 3 :    name_Value = " du: Navier Stokes equation ";
                        Value = min_u;
                        i_loc = i_u;
                        j_loc = j_u;
                        k_loc = k_u;
                        break;
        case 4 :    name_Value = " dv: Navier Stokes equation ";
                        Value = min_v;
                        i_loc = i_v;
                        j_loc = j_v;
                        k_loc = k_v;
                        break;
        case 5 :    name_Value = " dw: Navier Stokes equation ";
                        Value = min_w;
                        i_loc = i_w;
                        j_loc = j_w;
                        k_loc = k_w;
                        break;
        case 6 :    name_Value = " dt: energy transport equation ";
                        Value = min_t;
                        i_loc = i_t;
                        j_loc = j_t;
                        k_loc = k_t;
                        break;
        case 7 :    name_Value = " dh2o: h2o vap transport equation ";
                        Value = min_c;
                        i_loc = i_c;
                        j_loc = j_c;
                        k_loc = k_c;
                        break;
        case 8 :    name_Value = " dh2oc: h2o cl transport equation ";
                        Value = min_cloud;
                        i_loc = i_cloud;
                        j_loc = j_cloud;
                        k_loc = k_cloud;
                        break;
        case 9 :    name_Value = " dh2oi: h2o ic transport equation ";
                        Value = min_ice;
                        i_loc = i_ice;
                        j_loc = j_ice;
                        k_loc = k_ice;
                        break;
        case 10 :    name_Value = " dnh3: nh3 transport equation ";
                        Value = min_nh3;
                        i_loc = i_nh3;
                        j_loc = j_nh3;
                        k_loc = k_nh3;
                        break;
        case 11 :    name_Value = " dnh3: nh3 cl transport equation ";
                        Value = min_nh3_cloud;
                        i_loc = i_nh3_cloud;
                        j_loc = j_nh3_cloud;
                        k_loc = k_nh3_cloud;
                        break;
        case 12 :    name_Value = " dnh3i: nh3 i transport equation ";
                        Value = min_nh3_ice;
                        i_loc = i_nh3_ice;
                        j_loc = j_nh3_ice;
                        k_loc = k_nh3_ice;
                        break;
        case 13 :    name_Value = " dnh4sh: nh4sh transport equation ";
                        Value = min_nh4sh;
                        i_loc = i_nh4sh;
                        j_loc = j_nh4sh;
                        k_loc = k_nh4sh;
                        break;
        default :     cout << choice << "error in iterationPrintout_3D member function in class Accuracy" << endl;
    }
    i_loc_level = i_loc*int(L_atm)/(im-1)*1.e-3;
    if(j_loc <= 90){
        j_loc_deg = 90 - j_loc;
        deg_lat = deg_north;
    }
    if(j_loc > 90){
        j_loc_deg = j_loc - 90;
        deg_lat = deg_south;
    }
    if(k_loc <= 180){
        k_loc_deg = k_loc;
        deg_lon = deg_east;
    }
    if(k_loc > 180){
        k_loc_deg = 360 - k_loc;
        deg_lon = deg_west;
    }
    cout << setiosflags(ios::left) << setw(36) << setfill('.') << name_Value << " = " << resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << Value << setw(5) << j_loc_deg << setw(3) << deg_lat << setw(4) << k_loc_deg << setw(3) << deg_lon << setw(6) << i_loc_level << setw(2) << level << endl;
    choice++;
    if(choice <= 13) goto preparation;
    cout << endl << endl;
    return;
}

