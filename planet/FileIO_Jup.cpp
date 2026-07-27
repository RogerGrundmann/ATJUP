/*
 * Atmosphere General Circulation Modell(ATJUP) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in aa spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the boundary and initial conditions for diverse variables
*/
#include "cJupiterModel.h"
#include "Utils.h"

#include <cstdint>
#include <cstring>
#include <vector>

using namespace std;
using namespace JupiterUtils;

void cJupiterModel::writeData(){
    cout << endl << "      ATJUP: writeData" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    int i_radial = 20;
    paraview_vtk_radial(iter_n, i_radial);

    int j_longal = 112;
    paraview_vtk_longal(iter_n, j_longal);

    int k_zonal = 180;
    paraview_vtk_zonal(iter_n, k_zonal);


    if (paraview_panorama_vts_flag && panorama_cnt == panorama_print) {
        paraview_panorama_vts(iter_n);
//        paraview_sphere_vts(iter_n);
    }

    JupiterPlotData();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for writeData\n", elapsed.count() * 1e-9);

    cout << "      ATJUP: writeData ended" << endl;

    return;
}
/*
*
*/
void cJupiterModel::writeResults(){
    cout << endl << "      ATJUP: writeResults" << endl;
//    double coeff_mmWS = r_mix/r_h2o;  // coeff_mmWS = 1.2041/0.0094 [kg/m³/kg/m³] = 128,0827 [/]
//    double coeff_lv = lv_h2o /(cp_h2o * t_0);  // coefficient for the specific latent Evaporation heat(Condensation heat), coeff_lv = 9.1069 in [/]
//    double coeff_ls = ls /(cp_h2o * t_0);  // coefficient for the specific latent Evaporation heat(Condensation heat), coeff_ls = 10.3091 in [/]
//    double a, e;

    auto begin = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            precipitable_water.y[j][k] = 0.;  // precipitable water
        }
    }
    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            Q_Latent.x[0][j][k] = c43 * Q_Latent.x[1][j][k] - c13 * Q_Latent.x[2][j][k];
            Q_Latent.x[im-1][j][k] = c43 * Q_Latent.x[im-2][j][k] - c13 * Q_Latent.x[im-3][j][k];
            Q_Sensible.x[0][j][k] = c43 * Q_Sensible.x[1][j][k] - c13 * Q_Sensible.x[2][j][k];
            Q_Sensible.x[im-1][j][k] = c43 * Q_Sensible.x[im-2][j][k] - c13 * Q_Sensible.x[im-3][j][k];
            BuoyancyForce.x[0][j][k] = c43 * BuoyancyForce.x[1][j][k] - c13 * BuoyancyForce.x[2][j][k];
            BuoyancyForce.x[im-1][j][k] = c43 * BuoyancyForce.x[im-2][j][k] - c13 * BuoyancyForce.x[im-3][j][k];
/*
            nh3.x[0][j][k] = c43 * nh3.x[1][j][k] - c13 * nh3.x[2][j][k];
            nh3.x[im-1][j][k] = c43 * nh3.x[im-2][j][k] - c13 * nh3.x[im-3][j][k];
            nh3_cloud.x[0][j][k] = c43 * nh3_cloud.x[1][j][k] - c13 * nh3_cloud.x[2][j][k];
            nh3_cloud.x[im-1][j][k] = c43 * nh3_cloud.x[im-2][j][k] - c13 * nh3_cloud.x[im-3][j][k];
            nh3_ice.x[0][j][k] = c43 * nh3_ice.x[1][j][k] - c13 * nh3_ice.x[2][j][k];
            nh3_ice.x[im-1][j][k] = c43 * nh3_ice.x[im-2][j][k] - c13 * nh3_ice.x[im-3][j][k];
*/
        }
    }
    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            t.x[i][0][k] = c43 * t.x[i][1][k] - c13 * t.x[i][2][k];
            t.x[i][jm-1][k] = c43 * t.x[i][jm-2][k] - c13 * t.x[i][jm-3][k];
            Q_Latent.x[i][0][k] = c43 * Q_Latent.x[i][1][k] - c13 * Q_Latent.x[i][2][k];
            Q_Latent.x[i][jm-1][k] = c43 * Q_Latent.x[i][jm-2][k] - c13 * Q_Latent.x[i][jm-3][k];
            Q_Sensible.x[i][0][k] = c43 * Q_Sensible.x[i][1][k] - c13 * Q_Sensible.x[i][2][k];
            Q_Sensible.x[i][jm-1][k] = c43 * Q_Sensible.x[i][jm-2][k] - c13 * Q_Sensible.x[i][jm-3][k];
            BuoyancyForce.x[i][0][k] = c43 * BuoyancyForce.x[i][1][k] - c13 * BuoyancyForce.x[i][2][k];
            BuoyancyForce.x[i][jm-1][k] = c43 * BuoyancyForce.x[i][jm-2][k] - c13 * BuoyancyForce.x[i][jm-3][k];
/*
            nh3.x[i][0][k] = c43 * nh3.x[i][1][k] - c13 * nh3.x[i][2][k];
            nh3.x[i][jm-1][k] = c43 * nh3.x[i][jm-2][k] - c13 * nh3.x[i][jm-3][k];
            nh3_cloud.x[i][0][k] = c43 * nh3_cloud.x[i][1][k] - c13 * nh3_cloud.x[i][2][k];
            nh3_cloud.x[i][jm-1][k] = c43 * nh3_cloud.x[i][jm-2][k] - c13 * nh3_cloud.x[i][jm-3][k];
            nh3_ice.x[i][0][k] = c43 * nh3_ice.x[i][1][k] - c13 * nh3_ice.x[i][2][k];
            nh3_ice.x[i][jm-1][k] = c43 * nh3_ice.x[i][jm-2][k] - c13 * nh3_ice.x[i][jm-3][k];
            h2o.x[i][0][k] = c43 * h2o.x[i][1][k] - c13 * h2o.x[i][2][k];
            h2o.x[i][jm-1][k] = c43 * h2o.x[i][jm-2][k] - c13 * h2o.x[i][jm-3][k];
            h2o_cloud.x[i][0][k] = c43 * h2o_cloud.x[i][1][k] - c13 * h2o_cloud.x[i][2][k];
            h2o_cloud.x[i][jm-1][k] = c43 * h2o_cloud.x[i][jm-2][k] - c13 * h2o_cloud.x[i][jm-3][k];
            h2o_ice.x[i][0][k] = c43 * h2o_ice.x[i][1][k] - c13 * h2o_ice.x[i][2][k];
            h2o_ice.x[i][jm-1][k] = c43 * h2o_ice.x[i][jm-2][k] - c13 * h2o_ice.x[i][jm-3][k];
*/
/*
            P_rain.x[i][0][k] = c43 * P_rain.x[i][1][k] - c13 * P_rain.x[i][2][k];
            P_rain.x[i][jm-1][k] = c43 * P_rain.x[i][jm-2][k] - c13 * P_rain.x[i][jm-3][k];
            P_snow.x[i][0][k] = c43 * P_snow.x[i][1][k] - c13 * P_snow.x[i][2][k];
            P_snow.x[i][jm-1][k] = c43 * P_snow.x[i][jm-2][k] - c13 * P_snow.x[i][jm-3][k];
*/
        }
    }
    #pragma omp parallel for
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            t.x[i][j][0] = c43 * t.x[i][j][1] - c13 * t.x[i][j][2];
            t.x[i][j][km-1] = c43 * t.x[i][j][km-2] - c13 * t.x[i][j][km-3];
            t.x[i][j][0] = t.x[i][j][km-1] =(t.x[i][j][0] + t.x[i][j][km-1])/ 2.;
            Q_Latent.x[i][j][0] = c43 * Q_Latent.x[i][j][1] - c13 * Q_Latent.x[i][j][2];
            Q_Latent.x[i][j][km-1] = c43 * Q_Latent.x[i][j][km-2] - c13 * Q_Latent.x[i][j][km-3];
            Q_Latent.x[i][j][0] = Q_Latent.x[i][j][km-1] =(Q_Latent.x[i][j][0] + Q_Latent.x[i][j][km-1])/ 2.;
            Q_Sensible.x[i][j][0] = c43 * Q_Sensible.x[i][j][1] - c13 * Q_Sensible.x[i][j][2];
            Q_Sensible.x[i][j][km-1] = c43 * Q_Sensible.x[i][j][km-2] - c13 * Q_Sensible.x[i][j][km-3];
            Q_Sensible.x[i][j][0] = Q_Sensible.x[i][j][km-1] =(Q_Sensible.x[i][j][0] + Q_Sensible.x[i][j][km-1])/ 2.;
            BuoyancyForce.x[i][j][0] = c43 * BuoyancyForce.x[i][j][1] - c13 * BuoyancyForce.x[i][j][2];
            BuoyancyForce.x[i][j][km-1] = c43 * BuoyancyForce.x[i][j][km-2] - c13 * BuoyancyForce.x[i][j][km-3];
            BuoyancyForce.x[i][j][0] = BuoyancyForce.x[i][j][km-1] =(BuoyancyForce.x[i][j][0] + BuoyancyForce.x[i][j][km-1])/ 2.;
/*
            nh3.x[i][j][0] = c43 * nh3.x[i][j][1] - c13 * nh3.x[i][j][2];
            nh3.x[i][j][km-1] = c43 * nh3.x[i][j][km-2] - c13 * nh3.x[i][j][km-3];
            nh3.x[i][j][0] = nh3.x[i][j][km-1] =(nh3.x[i][j][0] + nh3.x[i][j][km-1])/ 2.;
            nh3_cloud.x[i][j][0] = c43 * nh3_cloud.x[i][j][1] - c13 * nh3_cloud.x[i][j][2];
            nh3_cloud.x[i][j][km-1] = c43 * nh3_cloud.x[i][j][km-2] - c13 * nh3_cloud.x[i][j][km-3];
            nh3_cloud.x[i][j][0] = nh3_cloud.x[i][j][km-1] =(nh3_cloud.x[i][j][0] + nh3_cloud.x[i][j][km-1])/ 2.;
            nh3_ice.x[i][j][0] = c43 * nh3_ice.x[i][j][1] - c13 * nh3_ice.x[i][j][2];
            nh3_ice.x[i][j][km-1] = c43 * nh3_ice.x[i][j][km-2] - c13 * nh3_ice.x[i][j][km-3];
            nh3_ice.x[i][j][0] = nh3_ice.x[i][j][km-1] =(nh3_ice.x[i][j][0] + nh3_ice.x[i][j][km-1])/ 2.;
            h2o_cloud.x[i][j][0] = c43 * h2o_cloud.x[i][j][1] - c13 * h2o_cloud.x[i][j][2];
            h2o_cloud.x[i][j][km-1] = c43 * h2o_cloud.x[i][j][km-2] - c13 * h2o_cloud.x[i][j][km-3];
            h2o_cloud.x[i][j][0] = h2o_cloud.x[i][j][km-1] =(h2o_cloud.x[i][j][0] + h2o_cloud.x[i][j][km-1])/ 2.;
            h2o_ice.x[i][j][0] = c43 * h2o_ice.x[i][j][1] - c13 * h2o_ice.x[i][j][2];
            h2o_ice.x[i][j][km-1] = c43 * h2o_ice.x[i][j][km-2] - c13 * h2o_ice.x[i][j][km-3];
            h2o_ice.x[i][j][0] = h2o_ice.x[i][j][km-1] =(h2o_ice.x[i][j][0] + h2o_ice.x[i][j][km-1])/ 2.;
*/
/*
            P_rain.x[i][j][0] = c43 * P_rain.x[i][j][1] - c13 * P_rain.x[i][j][2];
            P_rain.x[i][j][km-1] = c43 * P_rain.x[i][j][km-2] - c13 * P_rain.x[i][j][km-3];
            P_rain.x[i][j][0] = P_rain.x[i][j][km-1] =(P_rain.x[i][j][0] + P_rain.x[i][j][km-1])/ 2.;
            P_snow.x[i][j][0] = c43 * P_snow.x[i][j][1] - c13 * P_snow.x[i][j][2];
            P_snow.x[i][j][km-1] = c43 * P_snow.x[i][j][km-2] - c13 * P_snow.x[i][j][km-3];
            P_snow.x[i][j][0] = P_snow.x[i][j][km-1] =(P_snow.x[i][j][0] + P_snow.x[i][j][km-1])/ 2.;
*/
        }
    }
/*
    precipitable_water.y[0][0] = 0.;
    precipitable_water.y[0][0] = 0.;
    nh3_total.y[0][0] = 0.;
    nh3_cloud_total.y[0][0] = 0.;
    nh3_ice_total.y[0][0] = 0.;
    double precipitablewater_average = 0.;
    double precipitation_average = 0.;
    double h2_vegetation_average = 0.;
    double he_vegetation_average = 0.;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            nh3_total.y[j][k] = nh3.x[0][j][k];
            nh3_cloud_total.y[j][k] = nh3_cloud.x[0][j][k];
            nh3_ice_total.y[j][k] = nh3_ice.x[0][j][k];
            for(int i = 0; i < im; i++){
                e = h2o.x[i][j][k] * p_stat.x[i][j][k]/ep_h2o;  // water vapour pressure in hPa
                a = 216.6 * e /(t.x[i][j][k] * t_0);  // absolute humidity in kg/m3
                precipitable_water.y[j][k] += a * L_atm /(double)(im - 1);  // kg/m³ * m
            }
        }
    }
    double coeff_prec = 86400.;  // dimensions see below
// surface values of precipitation and precipitable water
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            // All-species surface total (H2O + NH3 + NH4SH), computed by PrecipitationJup at the
            // base of each column. This used to be H2O rain+snow only, read at i=0 and clipped at
            // 25 mm/d — an Earth-rainfall convention that ignored graupel, the whole NH3 deck and
            // the NH4SH crystals. The clip is dropped: the flux is already bounded by
            // PrecipitationJup's P_max_flux, and Jupiter's energy budget puts the physical values
            // near 1 mm/d, so a 25 mm/d ceiling only ever hid a runaway.
            Precipitation.y[j][k] = coeff_prec * precip_srf_total.y[j][k];
            // 60 s * 60 min * 24 h = 86400 s == 1 d
            // precip_srf_total in kg/ ( m² * s ) = mm/s
            // Precipitation in 86400. * kg/ ( m² * d ) = 86400 mm/d
            // kg/ ( m² * s ) == mm/s ( Kraus, p. 94 )
            if(Precipitation.y[j][k] <= 0)  Precipitation.y[j][k] = 0.;
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            precipitablewater_average += precipitable_water.y[j][k];
            precipitation_average += Precipitation.y[j][k];
            h2_vegetation_average += nh3_total.y[j][k];
            he_vegetation_average += nh3_cloud_total.y[j][k];
        }
    }
    h2_vegetation_average = h2_vegetation_average /(double)((jm-1)*(km-1));
    he_vegetation_average = he_vegetation_average /(double)((jm-1)*(km-1));
    precipitablewater_average = precipitablewater_average /(double)((jm-1)*(km-1));
    precipitation_average = 365. * precipitation_average /(double)((jm-1)*(km-1));
    cout.precision(2);
    string level = "m";
    string deg_north = "°N";
    string deg_south = "°S";
    string deg_west = "°W";
    string deg_east = "°E";
    string name_Value_7 = " precipitable water average ";
    string name_Value_8 = " precipitation average per year ";
    string name_Value_9 = " precipitation average per day ";
    string name_Value_22 = " H2_average ";
    string name_Value_25 = " He_average ";
    string name_unit_mmd = " mm/d";
    string name_unit_mm = " mm";
    string name_unit_mma = " mm/a";
    string name_unit_ppm = " kg/kg";
    cout << endl;
    double Value_7 = precipitablewater_average;
    double Value_8 = precipitation_average;
    cout << setw(6)<< setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_7 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_7 << setw(6)<< name_unit_mm << "   " << setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_8 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_8 << setw(6)<< name_unit_mma << "   " << setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_9 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_8/365. << setw(6)<< name_unit_mmd << endl;
    double Value_9 = h2_vegetation_average;
    double Value_25 = he_vegetation_average;
    cout << setw(6)<< setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_22 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_9 << setw(6)<< name_unit_ppm << "   " << setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_25 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_25 << setw(6)<< name_unit_ppm << endl << endl << endl;
    return;
*/

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for writeResults\n", elapsed.count() * 1e-9);

    cout << "      ATJUP: writeResults ended" << endl;

}
/*
*
*/


// ==================== FULL-3D-STATE CHECKPOINT / RESTART ====================
// Binary dump of the prognostic fields so a run can resume at a chosen iteration instead of
// re-spinning the circulation from scratch. The ATJUP counterpart of ATOM's
// cAtmosphereModel::save_state / load_state (FileIO_Atm.cpp), with the same file layout:
// a 5-int header followed by the arrays, each written as im*jm contiguous rows of km doubles.
//
// Only the genuinely PROGNOSTIC arrays are stored. Everything else — the reaction rates,
// diffusive and thermal mass fluxes, forces, radiation, precipitation fluxes and the latent
// and sensible heat fields — is recomputed from these at the top of every even iteration, so
// storing them would only add bulk and a way for the file to disagree with itself.
//
// p_stat is included even though it is quasi-static: it is what the buoyancy term and the
// whole saturation chain read, and it must match the temperature field it was built with.
std::vector<Array*> cJupiterModel::restart_arrays(){
    return { &t,   &u,   &v,   &w,
             &tn,  &un,  &vn,  &wn,
             &h2o,  &h2o_cloud,  &h2o_ice,
             &h2on, &h2o_cloudn, &h2o_icen,
             &h2s,  &h2sn,
             &nh3,  &nh3_cloud,  &nh3_ice,
             &nh3n, &nh3_cloudn, &nh3_icen,
             &ch4,  &ch4_cloud,  &ch4_ice,
             &ch4n, &ch4_cloudn, &ch4_icen,
             &nh4sh, &nh4shn,
             &p_dyn, &p_dynn, &p_stat,
             &tke, &dis, &tken, &disn, &nue };
}

// Names parallel to restart_arrays(), so the NaN report can say WHICH field went bad.
static const char* const restart_array_names[] = {
    "t","u","v","w",
    "tn","un","vn","wn",
    "h2o","h2o_cloud","h2o_ice",
    "h2on","h2o_cloudn","h2o_icen",
    "h2s","h2sn",
    "nh3","nh3_cloud","nh3_ice",
    "nh3n","nh3_cloudn","nh3_icen",
    "ch4","ch4_cloud","ch4_ice",
    "ch4n","ch4_cloudn","ch4_icen",
    "nh4sh","nh4shn",
    "p_dyn","p_dynn","p_stat",
    "tke","dis","tken","disn","nue" };

// Per-iteration NaN watch (ATJUP_NANCHECK=1). Scans the same fields the checkpoint stores
// and reports the FIRST non-finite cell with the field name and its (i,j,k), so a blow-up
// is caught at the iteration it starts rather than being discovered hundreds of iterations
// later in the output. Returns true while the state is clean.
bool cJupiterModel::nan_watch(int iter){
    std::vector<Array*> arrs = restart_arrays();
    const int n_names = (int)(sizeof(restart_array_names)/sizeof(restart_array_names[0]));

    // A CENSUS, not just the first cell in scan order. The scan runs i outermost, so a NaN
    // born in the interior gets reported at the i=0 boundary cell that merely inherited it
    // through the bcRadius extrapolation. Counting per field and recording the index extent
    // shows at a glance whether this is one interior cell, a whole boundary plane, or already
    // everywhere — which is the difference between a local physics problem and a global one.
    long total = 0;
    bool any = false;
    for(size_t a = 0; a < arrs.size(); a++){
        long n = 0;
        int i_lo = im, i_hi = -1, j_lo = jm, j_hi = -1, k_lo = km, k_hi = -1;
        int fi = -1, fj = -1, fk = -1;
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                for(int k = 0; k < km; k++){
                    std::uint64_t bits;
                    std::memcpy(&bits, &arrs[a]->x[i][j][k], sizeof(bits));
                    if((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL){
                        n++;
                        if(i < i_lo) i_lo = i;
                        if(i > i_hi) i_hi = i;
                        if(j < j_lo) j_lo = j;
                        if(j > j_hi) j_hi = j;
                        if(k < k_lo) k_lo = k;
                        if(k > k_hi) k_hi = k;
                        if(fi < 0){ fi = i; fj = j; fk = k; }
                    }
                }
        if(n > 0){
            if(!any){
                printf("      ATJUP: ===== NAN WATCH: state went non-finite at iteration %d =====\n",
                       iter);
                any = true;
            }
            printf("        %-12s %8ld cells   i[%d..%d] j[%d..%d] k[%d..%d]   first (%d,%d,%d)"
                   "  SeaMount=%g t=%g p_stat=%g\n",
                   (a < (size_t)n_names ? restart_array_names[a] : "?"), n,
                   i_lo, i_hi, j_lo, j_hi, k_lo, k_hi, fi, fj, fk,
                   SeaMount.x[fi][fj][fk], t.x[fi][fj][fk], p_stat.x[fi][fj][fk]);
            total += n;
        }
    }
    if(any) printf("        total %ld non-finite cells\n", total);
    return !any;
}

// Per-radial-level momentum census (ATJUP_WPROFILE=<stride>, 0 = off).
//
// printMinMax only reports ONE global extremum per field, which cannot distinguish a
// boundary artifact from a domain-wide momentum source: both look like "max|w| is growing".
// This prints, for every level i, the level maximum and the area-weighted level mean of each
// velocity component, so the VERTICAL STRUCTURE of the growth is visible:
//   growth confined to i = 1..3      -> the radial boundary treatment is making it,
//   growth at every level together   -> the momentum budget has no sink and it is physical.
// The level mean is the discriminating number: an extrapolation artifact raises the extremum
// near the wall without moving the layer's mean momentum, a missing sink raises both.
//
// One line per level, prefixed WPROF and space-separated, so a run log can be reduced with
// grep/awk without parsing the surrounding report. Velocities are printed in m/s.
void cJupiterModel::momentum_profile(int iter){
    printf("      ATJUP: ===== MOMENTUM PROFILE at iteration %d =====\n", iter);
    printf("      WPROF iter    i   z[km]     max|u|     max|v|     max|w|"
           "      <u>        <v>        <w>      rms(w)   jmax kmax\n");

    for(int i = 0; i < im; i++){
        double mu = 0.0, mv = 0.0, mw = 0.0;
        double su = 0.0, sv = 0.0, sw = 0.0, sww = 0.0, wsum = 0.0;
        int bj = -1, bk = -1;

        for(int j = 0; j < jm; j++){
            const double wgt = sin(the.z[j]);          // spherical area weight, as in computeBuoyancyRefLevel
            for(int k = 0; k < km; k++){
                if(SeaMount.x[i][j][k] == 1.0) continue;      // solid cell: not part of the fluid budget
                const double uu = u.x[i][j][k], vv = v.x[i][j][k], ww = w.x[i][j][k];
                if(!std::isfinite(uu) || !std::isfinite(vv) || !std::isfinite(ww)) continue;
                if(fabs(uu) > mu) mu = fabs(uu);
                if(fabs(vv) > mv) mv = fabs(vv);
                if(fabs(ww) > mw){ mw = fabs(ww); bj = j; bk = k; }
                su   += wgt * uu;
                sv   += wgt * vv;
                sw   += wgt * ww;
                sww  += wgt * ww * ww;
                wsum += wgt;
            }
        }
        const double n = (wsum > 0.0) ? wsum : 1.0;
        printf("      WPROF %4d %4d %7.2f %10.3f %10.3f %10.3f %10.4f %10.4f %10.4f %10.4f %5d %5d\n",
               iter, i, get_layer_height(i) , mu * u_0, mv * u_0, mw * u_0,
               (su / n) * u_0, (sv / n) * u_0, (sw / n) * u_0,
               sqrt(sww / n) * u_0, bj, bk);
    }
}

// Floor the condensable species at zero, and keep the books on what that costs.
//
// A negative concentration has no meaning, and the fields do go slightly negative: the central
// differences of the transport terms undershoot wherever a species has a sharp edge, and the
// (4/3,-1/3) extrapolation at the radial boundary planes undershoots a field that is already
// essentially zero there. Measured at iteration 200 the clipped amount is small enough that a
// plain floor is the right answer rather than a mass-conserving filler:
//
//   nh4sh      0.017 % of the positive mass, 125418 cells
//   nh3_cloud  0.0013 %      nh3_ice 0.0009 %      h2o_ice 0.0001 %
//   h2o, h2s, ch4 and the ch4 condensates: zero or a handful of cells at 1e-11
//
// and 96 % of those nh4sh cells sit on the two radial boundary planes i=0 and i=40, where the
// field is 1e-16 and smaller — they are extrapolation noise, not transport undershoot. The
// genuine interior ones cluster at i=10, the widest flank of the obstacle.
//
// A floor is in principle a mass SOURCE, so the clipped amount is accumulated per field and
// reported next to printMinMax rather than left invisible. Read that report as GROSS clipping,
// not as net mass gained: h2o_cloud accumulates 19.6 % of its own mass in 50 iterations, which
// looks alarming and is not, because SaturationAdjustmentJup re-partitions vapour and condensate
// on the next pass and gives it straight back. Verified by running the same 200 iterations with
// and without the floor and summing h2o + h2o_cloud + h2o_ice:
//
//   floor on   total water  3.3293e4 -> 3.2436e4   (-2.573 % over 100 iterations)
//   floor off  total water  3.3279e4 -> 3.2432e4   (-2.545 %)
//
// i.e. the two budgets agree to 0.03 percentage points, while h2o_cloud's minimum goes from
// -5.07e-4 to exactly 0. The floor buys a clean field for the price of nothing measurable.
// What the counter is FOR is the day that stops being true. ATJUP_NO_CLAMP=1 turns it off.
void cJupiterModel::clampNegativeSpecies(){
    static const int off = [](){ const char* e = getenv("ATJUP_NO_CLAMP"); return e ? atoi(e) : 0; }();
    if(off) return;

    Array* fields[] = {
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s,
        &nh3, &nh3_cloud, &nh3_ice,
        &ch4, &ch4_cloud, &ch4_ice,
        &nh4sh };
    const int nf = (int)(sizeof(fields) / sizeof(fields[0]));

    if((int)clamp_added.size() != nf){
        clamp_added.assign(nf, 0.0);
        clamp_cells.assign(nf, 0);
    }

    for(int f = 0; f < nf; f++){
        Array& F = *fields[f];
        double added = 0.0;
        long   cells = 0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:added,cells)
        for(int i = 0; i < im; i++){
            for(int j = 0; j < jm; j++){
                for(int k = 0; k < km; k++){
                    const double v = F.x[i][j][k];
                    // Written as !(v >= 0.0) so a NaN is caught here too rather than carried on.
                    if(!(v >= 0.0)){
                        if(std::isfinite(v)){ added -= v; cells++; }
                        F.x[i][j][k] = 0.0;
                    }
                }
            }
        }
        clamp_added[f] += added;
        clamp_cells[f] += cells;
    }
}

// Companion report, called from printMinMax so it shares the checkpoint cadence.
void cJupiterModel::reportClampBudget(){
    static const char* const names[] = {
        "h2o","h2o_cloud","h2o_ice","h2s","nh3","nh3_cloud","nh3_ice",
        "ch4","ch4_cloud","ch4_ice","nh4sh" };
    const int nf = (int)(sizeof(names)/sizeof(names[0]));
    if((int)clamp_added.size() != nf) return;

    Array* fields[] = {
        &h2o, &h2o_cloud, &h2o_ice, &h2s, &nh3, &nh3_cloud, &nh3_ice,
        &ch4, &ch4_cloud, &ch4_ice, &nh4sh };

    bool any = false;
    for(int f = 0; f < nf; f++) if(clamp_cells[f] > 0) any = true;
    if(!any) return;

    printf("\n      ATJUP: negative-value clamp, cumulative since start\n");
    for(int f = 0; f < nf; f++){
        if(clamp_cells[f] == 0) continue;
        double pos = 0.0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:pos)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                for(int k = 0; k < km; k++)
                    if(fields[f]->x[i][j][k] > 0.0) pos += fields[f]->x[i][j][k];
        printf("        %-12s gross %.4e over %10ld clippings = %8.4f %% of the current"
               " field mass (gross, not net — see the note in FileIO_Jup.cpp)\n",
               names[f], clamp_added[f], clamp_cells[f],
               (pos > 0.0) ? 100.0 * clamp_added[f] / pos : 0.0);
    }
}

void cJupiterModel::save_state(int iter){
    const string fn = output_path + "/jup_restart_" + std::to_string(iter) + ".bin";
    std::ofstream f(fn, std::ios::binary);
    if(!f){
        cout << "      ATJUP: save_state FAILED to open " << fn << endl;
        return;
    }
    // Header: magic, grid dimensions, and the iteration this state belongs to. The grid is
    // checked on load so a restart written at a different resolution is rejected rather than
    // read as garbage.
    const int32_t hdr[5] = { 0x4A555031 /*"JUP1"*/, im, jm, km, iter };
    f.write(reinterpret_cast<const char*>(hdr), sizeof(hdr));

    std::vector<Array*> arrs = restart_arrays();
    for(size_t a = 0; a < arrs.size(); a++)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                f.write(reinterpret_cast<const char*>(arrs[a]->x[i][j]), km * sizeof(double));

    if(!f){
        cout << "      ATJUP: save_state FAILED while writing " << fn
             << " (disk full?)" << endl;
        return;
    }
    const double mb = (double)(sizeof(hdr) + arrs.size() * (size_t)im * jm * km * sizeof(double))
                    / (1024.0 * 1024.0);
    printf("      ATJUP: save_state wrote %zu arrays (%.1f MB) to %s\n",
           arrs.size(), mb, fn.c_str());
}

bool cJupiterModel::load_state(int iter){
    const string fn = output_path + "/jup_restart_" + std::to_string(iter) + ".bin";
    std::ifstream f(fn, std::ios::binary);
    if(!f){
        cout << "      ATJUP: load_state: no file " << fn
             << " - running from scratch" << endl;
        return false;
    }
    int32_t hdr[5];
    f.read(reinterpret_cast<char*>(hdr), sizeof(hdr));
    if(!f || hdr[0] != 0x4A555031 || hdr[1] != im || hdr[2] != jm || hdr[3] != km){
        cout << "      ATJUP: load_state: bad header / grid mismatch in " << fn
             << " - running from scratch" << endl;
        return false;
    }

    std::vector<Array*> arrs = restart_arrays();
    for(size_t a = 0; a < arrs.size(); a++)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++){
                f.read(reinterpret_cast<char*>(arrs[a]->x[i][j]), km * sizeof(double));
                if(!f){
                    cout << "      ATJUP: load_state: truncated file " << fn
                         << " - running from scratch" << endl;
                    return false;
                }
            }

    cout << "      ATJUP: load_state restored " << arrs.size() << " arrays from "
         << fn << " (resuming after iteration " << hdr[4] << ")" << endl;
    return true;
}

// True when every serialized prognostic field is finite everywhere. Guards the periodic
// checkpoint: a diverged state must never overwrite a good restart point, because the whole
// value of the file is that you can resume from it. Uses the IEEE-754 exponent bits rather
// than std::isfinite for the reason given in RungeKutta_Jup_Turb.cpp.
bool cJupiterModel::restart_state_is_clean(){
    std::vector<Array*> arrs = restart_arrays();
    bool clean = true;
    for(size_t a = 0; a < arrs.size() && clean; a++)
        for(int i = 0; i < im && clean; i++)
            for(int j = 0; j < jm && clean; j++)
                for(int k = 0; k < km; k++){
                    std::uint64_t bits;
                    std::memcpy(&bits, &arrs[a]->x[i][j][k], sizeof(bits));
                    if((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL){
                        clean = false;
                        break;
                    }
                }
    return clean;
}
/*
*
*/
