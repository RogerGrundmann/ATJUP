/*
 * Jupiter General Circulation Modell(ATJUP)applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write sequel, transfer and paraview files
*/

#include "cJupiterModel.h"
#include "ParaViewWriter.h"

using namespace std;
//using namespace JupiterUtils;

// The five dumpers moved to the SHARED ParaViewWriter.h, as namespace ParaViewIO.
// They were already byte-identical in both models; see that header.
/*
 * 
*/
void cJupiterModel::paraview_panorama_vts(int n){
    using namespace ParaViewIO;
    double x, y, z, dx, dy, dz;
    double r_mix_plus = r_mix * 1e6;
    string Jupiter_panorama_vts_File_Name = output_path + "/Jupiter_panorama_" 
        + std::to_string(n) + ".vts";
    ofstream Jupiter_panorama_vts_File;
    Jupiter_panorama_vts_File.precision(4);
    Jupiter_panorama_vts_File.setf(ios::fixed);
    Jupiter_panorama_vts_File.open(Jupiter_panorama_vts_File_Name);
    if(!Jupiter_panorama_vts_File.is_open()){
        cerr << "ERROR: could not open shpere_vts file " << __FILE__ 
            << " at line " << __LINE__ << "\n";
        abort();
    }
    Jupiter_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Jupiter_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Jupiter_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Jupiter_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Jupiter_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Temperature PressureDynamic PressureStatic CH4 CH4Cloud CH4Ice NH3 NH3Cloud NH3Ice H2O H2OCloud H2OIce Q_Latent Q_Sensible Q_rad_mW_m3 Radiation P_rain_mmd P_snow_mmd P_graupel_mmd P_nh3_rain_mmd P_nh3_snow_mmd P_nh3_graupel_mmd P_ch4_rain_mmd P_ch4_snow_mmd P_ch4_graupel_mmd P_nh4sh_mmd Q_precip_mW_m3 tke_m2s2 dis_nd nue_t_m2s prod_nd tke_source_nd dis_source_nd BuoyancyForce \">\n"  << endl;

    Jupiter_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_panorama_vts_File << u.x[i][j][k] << " " << v.x[i][j][k] << " " << w.x[i][j][k] << endl;
            }
            Jupiter_panorama_vts_File <<  "\n"  << endl;
        }
        Jupiter_panorama_vts_File <<  "\n"  << endl;
    }
    Jupiter_panorama_vts_File <<  "\n"  << endl;
    Jupiter_panorama_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_panorama_vts_File << t.x[i][j][k] * t_ref - 273.15 << endl;
            }
            Jupiter_panorama_vts_File <<  "\n"  << endl;
        }
        Jupiter_panorama_vts_File <<  "\n"  << endl;
    }
    Jupiter_panorama_vts_File <<  "\n"  << endl;
    Jupiter_panorama_vts_File <<  "    </DataArray>\n" << endl;
    dump_array("Seamount", SeaMount, 1.0, Jupiter_panorama_vts_File);
    dump_array("u-component", u, u_0, Jupiter_panorama_vts_File);
    dump_array("v-component", v, u_0, Jupiter_panorama_vts_File);
    dump_array("w-component", w, u_0, Jupiter_panorama_vts_File);

    // PressureDyn goes out in MILLIBAR. p_dyn is stored as the nondimensional kinematic
    // pressure (see p_dyn_to_bar in cJupiterModel.h); in bar it would be ~0.002, and these
    // writers print four decimals, so bar would throw away all but one digit. Same reason the
    // precipitation fluxes go out in mm/day and Q_rad in mW/m3.
    dump_array("PressureDyn", p_dyn, p_dyn_to_bar() * 1.0e3, Jupiter_panorama_vts_File);
//    dump_array("PressureStat", p_stat, 1.0, Jupiter_panorama_vts_File);

//    dump_array("CoriolisForce", CoriolisForce, 1.0, Jupiter_panorama_vts_File);
//    dump_array("CentrifugalForce", CentrifugalForce, 1e3, Jupiter_panorama_vts_File);
//    dump_array("BuoyancyForce", BuoyancyForce, 1.0, Jupiter_panorama_vts_File);
//    dump_array("PresGradForce", PresGradForce, 1.0, Jupiter_panorama_vts_File);

    dump_array("CH4", ch4, 1.0, Jupiter_panorama_vts_File);
    dump_array("CH4Cloud", ch4_cloud, 1.0, Jupiter_panorama_vts_File);
    dump_array("CH4Ice", ch4_ice, 1.0, Jupiter_panorama_vts_File);

    dump_array("H2O", h2o, 1.0, Jupiter_panorama_vts_File);
    dump_array("H2OCloud", h2o_cloud, 1.0, Jupiter_panorama_vts_File);
    dump_array("H2OIce", h2o_ice, 1.0, Jupiter_panorama_vts_File);

    dump_array("H2S", h2s, 1.0, Jupiter_panorama_vts_File);
//    dump_array("w_h2s", w_h2s, 1.0, Jupiter_panorama_vts_File);

    dump_array("NH3", nh3, 1.0, Jupiter_panorama_vts_File);
    dump_array("NH3Cloud", nh3_cloud, 1.0, Jupiter_panorama_vts_File);
    dump_array("NH3Ice", nh3_ice, 1.0, Jupiter_panorama_vts_File);
//    dump_array("w_nh3", w_nh3, 1.0, Jupiter_panorama_vts_File);
//    dump_array("j_nh3", j_nh3, 1.0, Jupiter_panorama_vts_File);
//    dump_array("jT_nh3", jT_nh3, 1.0, Jupiter_panorama_vts_File);

    dump_array("NH4SH", nh4sh, r_mix_plus, Jupiter_panorama_vts_File);
//    dump_array("w_nh4sh", w_nh4sh, 1.0, Jupiter_panorama_vts_File);

    dump_array("Q_Latent", Q_Latent, 1.0, Jupiter_panorama_vts_File);
//    dump_array("Q_Sensible", Q_Sensible, 1.0, Jupiter_panorama_vts_File);

    dump_array("Q_rad_mW_m3", Q_rad, 1.0e3, Jupiter_panorama_vts_File);
    dump_array("Radiation", radiation, 1.0, Jupiter_panorama_vts_File);

    // Precipitation fluxes from PrecipitationJup, written as mm/day (x86400 from kg/m2/s;
    // 1 kg/m2 == 1 mm depth). H2O and NH3 each carry a full rain/snow/graupel triple; NH4SH
    // settles as crystals only. NOTE these streams use precision(4)+ios::fixed, so the raw SI
    // values (~1e-6 kg/m2/s) would every one of them round to 0.0000 — hence the unit scaling,
    // and the _mmd / _mW_m3 suffixes so the ParaView legend states the units.
    dump_array("P_rain_mmd", P_rain, 86400.0, Jupiter_panorama_vts_File);
    dump_array("P_snow_mmd", P_snow, 86400.0, Jupiter_panorama_vts_File);
    dump_array("P_graupel_mmd", P_graupel, 86400.0, Jupiter_panorama_vts_File);
    dump_array("P_nh3_rain_mmd", P_nh3_rain, 86400.0, Jupiter_panorama_vts_File);
    dump_array("P_nh3_snow_mmd", P_nh3_snow, 86400.0, Jupiter_panorama_vts_File);
    dump_array("P_nh3_graupel_mmd", P_nh3_graupel, 86400.0, Jupiter_panorama_vts_File);
    dump_array("P_ch4_rain_mmd", P_ch4_rain, 86400.0, Jupiter_panorama_vts_File);
    dump_array("P_ch4_snow_mmd", P_ch4_snow, 86400.0, Jupiter_panorama_vts_File);
    dump_array("P_ch4_graupel_mmd", P_ch4_graupel, 86400.0, Jupiter_panorama_vts_File);
    dump_array("P_nh4sh_mmd", P_nh4sh, 86400.0, Jupiter_panorama_vts_File);
    dump_array("Q_precip_mW_m3", Q_precip, 1.0e3, Jupiter_panorama_vts_File);

    // Turbulence (TurbulenceJup: k-epsilon / k-omega / k-omega SST). Zero unless ATJUP_TURB is
    // set. These streams use precision(4)+ios::fixed, so anything below 5e-5 rounds to 0.0000:
    // nue* is ~7e-5 and would vanish, and k* ~7e-3 would keep only two digits. The two worst-hit
    // fields are therefore written in PHYSICAL units, with the unit in the name:
    //   tke_m2s2  = k*   * u_0^2            [m2/s2]  (~67 for k* = 6.7e-3)
    //   nue_t_m2s = nue* * u_0 * L_atm[m]   [m2/s]   (~994, i.e. at ATOM's 1000 cap)
    // dis stays DIMENSIONLESS: its conversion is model-dependent (eps* uses u_0^3/L_atm, omega*
    // uses u_0/L_atm) and the writer cannot know which model ran; its raw range (3.6e-3 for
    // k-eps, 19 for k-omega) resolves fine. prod and the two source terms are left dimensionless
    // for the same reason — their magnitudes (1e-2 .. 1e2) print without loss.
    dump_array("tke_m2s2", tke, u_0 * u_0, Jupiter_panorama_vts_File);
    dump_array("dis_nd", dis, 1.0, Jupiter_panorama_vts_File);
    dump_array("nue_t_m2s", nue, u_0 * L_atm * 1.0e3, Jupiter_panorama_vts_File);
    dump_array("prod_nd", prod, 1.0, Jupiter_panorama_vts_File);
    dump_array("tke_source_nd", tke_source, 1.0, Jupiter_panorama_vts_File);
    dump_array("dis_source_nd", dis_source, 1.0, Jupiter_panorama_vts_File);

    Jupiter_panorama_vts_File <<  "   </PointData>\n" << endl;
    Jupiter_panorama_vts_File <<  "   <Points>\n"  << endl;
    Jupiter_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.1;
    dz = 0.1;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(k == 0 || j == 0) x = 0.0;
                else x = x + dx;
                Jupiter_panorama_vts_File << x << " " << y << " " << z  << endl;
            }
            x = 0;
            y = y + dy;
            Jupiter_panorama_vts_File <<  "\n"  << endl;
        }
        y = 0.0;
        z = z + dz;
        Jupiter_panorama_vts_File <<  "\n"  << endl;
    }
    Jupiter_panorama_vts_File <<  "    </DataArray>\n"  << endl;
    Jupiter_panorama_vts_File <<  "   </Points>\n"  << endl;
    Jupiter_panorama_vts_File <<  "  </Piece>\n"  << endl;
    Jupiter_panorama_vts_File <<  " </StructuredGrid>\n"  << endl;
    Jupiter_panorama_vts_File <<  "</VTKFile>\n"  << endl;
    Jupiter_panorama_vts_File.close();
    cout << "   File:  " << "Jupiter_panorama_" 
        << n << ".vts" << "  has been written to Directory:  " 
        << output_path << endl;
    return;
}
/*
 * 
*/
void cJupiterModel::paraview_vtk_radial(int n, int i_radial){
    using namespace ParaViewIO;
    double r_mix_plus = r_mix * 1e6;
    ofstream Jupiter_vtk_radial_File = ParaViewWriter<cJupiterModel>(*this)
        .open_slice("radial", "Radial", i_radial, n, km, jm, 0.1, false);
    const double z = 0.0;   // out-of-plane component of the in-plane vector below
    dump_radial("Seamount", SeaMount, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("u-Component", u, u_0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("v-Component", v, u_0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("w-Component", w, u_0, i_radial, Jupiter_vtk_radial_File);
    Jupiter_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
    Jupiter_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Jupiter_vtk_radial_File << t.x[i_radial][j][k] * t_ref - 273.15 << endl;
        }
    }

    dump_radial("thermalmassflux", thermalmassflux, 1e-3, i_radial, Jupiter_vtk_radial_File);
    dump_radial("CH4", ch4, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("CH4Cloud", ch4_cloud, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("CH4Ice", ch4_ice, 1.0, i_radial, Jupiter_vtk_radial_File);

    dump_radial("H2O", h2o, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("H2OCloud", h2o_cloud, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("H2OIce", h2o_ice, 1.0, i_radial, Jupiter_vtk_radial_File);

    dump_radial("H2S", h2s, 1.0, i_radial, Jupiter_vtk_radial_File);
/*
    dump_radial("w_h2s", w_h2s, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("j_h2s", j_h2s, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("jT_h2s", jT_h2s, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("massflux_h2s", massflux_h2s, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("difflux_h2s", difflux_h2s, 1.0, i_radial, Jupiter_vtk_radial_File);
*/
    dump_radial("NH3", nh3, 1.0, i_radial, Jupiter_vtk_radial_File);
/*
    dump_radial("NH3Cloud", nh3_cloud, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("NH3Ice", nh3_ice, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("w_nh3", w_nh3, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("j_nh3", j_nh3, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("jT_nh3", jT_nh3, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("massflux_nh3", massflux_nh3, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("difflux_nh3", difflux_nh3, 1.0, i_radial, Jupiter_vtk_radial_File);
*/
    dump_radial("NH4SH", nh4sh, r_mix_plus, i_radial, Jupiter_vtk_radial_File);
/*
    dump_radial("w_nh4sh", w_nh4sh, r_mix_plus, i_radial, Jupiter_vtk_radial_File);
    dump_radial("massflux_nh4sh", massflux_nh4sh, r_mix_plus, i_radial, Jupiter_vtk_radial_File);
    dump_radial("j_nh4sh", j_nh4sh, r_mix_plus, i_radial, Jupiter_vtk_radial_File);
    dump_radial("jT_nh4sh", jT_nh4sh, r_mix_plus, i_radial, Jupiter_vtk_radial_File);
    dump_radial("difflux_nh4sh", difflux_nh4sh, r_mix_plus, i_radial, Jupiter_vtk_radial_File);
*/

    dump_radial("PressureDyn", p_dyn, p_dyn_to_bar() * 1.0e3, i_radial, Jupiter_vtk_radial_File);
    dump_radial("PressureStat", p_stat, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("rho_mix", rho_mix, 1.0, i_radial, Jupiter_vtk_radial_File);

    dump_radial("CoriolisForce", CoriolisForce, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("CentrifugalForce", CentrifugalForce, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("BuoyancyForce", BuoyancyForce, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("PresGradForce", PresGradForce, 1.0, i_radial, Jupiter_vtk_radial_File);

    dump_radial("Q_Latent", Q_Latent, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("Q_Sensible", Q_Sensible, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("Q_rad", Q_rad, 1.0e3, i_radial, Jupiter_vtk_radial_File);
    dump_radial("Radiation", radiation, 1.0, i_radial, Jupiter_vtk_radial_File);

    // Precipitation fluxes, written as mm/day (x86400 from kg/m2/s; 1 kg/m2 == 1 mm depth).
    // These streams use precision(4)+ios::fixed, so raw SI values ~1e-6 would all round to
    // 0.0000 — hence the unit scaling, and the _mmd suffix so the legend is unambiguous.
    // Q_precip/Q_rad are scaled to mW/m3 for the same reason. Zero unless ATJUP_PRECIP is set.
    dump_radial("P_rain", P_rain, 86400.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("P_snow", P_snow, 86400.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("P_graupel", P_graupel, 86400.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("P_nh3_rain", P_nh3_rain, 86400.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("P_nh3_snow", P_nh3_snow, 86400.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("P_nh3_graupel", P_nh3_graupel, 86400.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("P_ch4_rain", P_ch4_rain, 86400.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("P_ch4_snow", P_ch4_snow, 86400.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("P_ch4_graupel", P_ch4_graupel, 86400.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("P_nh4sh", P_nh4sh, 86400.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("Q_precip", Q_precip, 1.0e3, i_radial, Jupiter_vtk_radial_File);

    // Turbulence (TurbulenceJup: k-epsilon / k-omega / k-omega SST). Zero unless ATJUP_TURB is
    // set. These streams use precision(4)+ios::fixed, so anything below 5e-5 rounds to 0.0000:
    // nue* is ~7e-5 and would vanish, and k* ~7e-3 would keep only two digits. The two worst-hit
    // fields are therefore written in PHYSICAL units, with the unit in the name:
    //   tke_m2s2  = k*   * u_0^2            [m2/s2]  (~67 for k* = 6.7e-3)
    //   nue_t_m2s = nue* * u_0 * L_atm[m]   [m2/s]   (~994, i.e. at ATOM's 1000 cap)
    // dis stays DIMENSIONLESS: its conversion is model-dependent (eps* uses u_0^3/L_atm, omega*
    // uses u_0/L_atm) and the writer cannot know which model ran; its raw range (3.6e-3 for
    // k-eps, 19 for k-omega) resolves fine. prod and the two source terms are left dimensionless
    // for the same reason — their magnitudes (1e-2 .. 1e2) print without loss.
    dump_radial("tke", tke, u_0 * u_0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("disd", dis, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("nue_t", nue, u_0 * L_atm * 1.0e3, i_radial, Jupiter_vtk_radial_File);
    dump_radial("prod", prod, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("tke_source", tke_source, 1.0, i_radial, Jupiter_vtk_radial_File);
    dump_radial("dis_source", dis_source, 1.0, i_radial, Jupiter_vtk_radial_File);

    // All-species SURFACE precipitation map, in mm/day. Unlike the fields above — which are
    // sampled at this file's fixed altitude i_radial and so miss any deck that does not live
    // there — these are 2D and always taken at the base of each column, giving the total
    // condensate mass flux actually arriving at the bottom plus its per-species breakdown.
    dump_radial_2d("Precip_total", precip_srf_total, 86400.0, Jupiter_vtk_radial_File);
    dump_radial_2d("Precip_h2o",   precip_srf_h2o,   86400.0, Jupiter_vtk_radial_File);
    dump_radial_2d("Precip_nh3",   precip_srf_nh3,   86400.0, Jupiter_vtk_radial_File);
    dump_radial_2d("Precip_ch4",   precip_srf_ch4,   86400.0, Jupiter_vtk_radial_File);
    dump_radial_2d("Precip_nh4sh", precip_srf_nh4sh, 86400.0, Jupiter_vtk_radial_File);

    // Per-column friction velocity u_tau from TurbulenceJup::compute_vel_star.
    dump_radial_2d("vel_star_ms", vel_star, 1.0, Jupiter_vtk_radial_File);

    Jupiter_vtk_radial_File <<  "VECTORS v-w-Cell float " << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Jupiter_vtk_radial_File << v.x[i_radial][j][k] << " " << w.x[i_radial][j][k] << " " << z << endl;
        }
    }
    ParaViewWriter<cJupiterModel>(*this).close_slice(Jupiter_vtk_radial_File, "radial", i_radial, n);
    return;
}
/*
 * 
*/
void cJupiterModel::paraview_vtk_zonal(int n, int k_zonal){
    using namespace ParaViewIO;
    double r_mix_plus = r_mix * 1e6;
    ofstream Jupiter_vtk_zonal_File = ParaViewWriter<cJupiterModel>(*this)
        .open_slice("zonal", "Zonal", k_zonal, n, jm, im, 0.05, false);
    const double z = 0.0;   // out-of-plane component of the in-plane vector below
    dump_zonal("Seamount", SeaMount, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("u-Component", u, u_0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("v-Component", v, u_0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("w-Component", w, u_0, k_zonal, Jupiter_vtk_zonal_File);
    Jupiter_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
    Jupiter_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;

    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Jupiter_vtk_zonal_File << t.x[i][j][k_zonal] * t_ref - 273.15 << endl;
            aux.x[i][j][k_zonal] = get_layer_height(i);
        }
    }
    dump_zonal("thermalmassflux", thermalmassflux, 1e-3, k_zonal, Jupiter_vtk_zonal_File);

    dump_zonal("height", aux, 1.0, k_zonal, Jupiter_vtk_zonal_File);

    dump_zonal("CH4", ch4, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("CH4Cloud", ch4_cloud, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("CH4Ice", ch4_ice, 1.0, k_zonal, Jupiter_vtk_zonal_File);

    dump_zonal("H2O", h2o, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("H2OCloud", h2o_cloud, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("H2OIce", h2o_ice, 1.0, k_zonal, Jupiter_vtk_zonal_File);

    dump_zonal("H2S", h2s, 1.0, k_zonal, Jupiter_vtk_zonal_File);
/*
    dump_zonal("w_h2s", w_h2s, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("j_h2s", j_h2s, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("jT_h2s", jT_h2s, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("massflux_h2s", massflux_h2s, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("difflux_h2s", difflux_h2s, 1.0, k_zonal, Jupiter_vtk_zonal_File);
*/
    dump_zonal("NH3", nh3, 1.0, k_zonal, Jupiter_vtk_zonal_File);
/*
    dump_zonal("NH3Cloud", nh3_cloud, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("NH3Ice", nh3_ice, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("w_nh3", w_nh3, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("j_nh3", j_nh3, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("jT_nh3", jT_nh3, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("massflux_nh3", massflux_nh3, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("difflux_nh3", difflux_nh3, 1.0, k_zonal, Jupiter_vtk_zonal_File);
*/
    dump_zonal("NH4SH", nh4sh, r_mix_plus, k_zonal, Jupiter_vtk_zonal_File);
/*
    dump_zonal("w_nh4sh", w_nh4sh, r_mix_plus, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("massflux_nh4sh", massflux_nh4sh, r_mix_plus, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("j_nh4sh", j_nh4sh, r_mix_plus, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("jT_nh4sh", jT_nh4sh, r_mix_plus, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("difflux_nh4sh", difflux_nh4sh, r_mix_plus, k_zonal, Jupiter_vtk_zonal_File);
*/
    dump_zonal("PressureDyn", p_dyn, p_dyn_to_bar() * 1.0e3, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("PressureStat", p_stat, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("rho_mix", rho_mix, 1.0, k_zonal, Jupiter_vtk_zonal_File);

    dump_zonal("CoriolisForce", CoriolisForce, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("CentrifugalForce", CentrifugalForce, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("BuoyancyForce", BuoyancyForce, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("PresGradForce", PresGradForce, 1.0, k_zonal, Jupiter_vtk_zonal_File);

    dump_zonal("Q_Latent", Q_Latent, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("Q_Sensible", Q_Sensible, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("Q_rad", Q_rad, 1.0e3, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("Radiation", radiation, 1.0, k_zonal, Jupiter_vtk_zonal_File);

    // Precipitation fluxes, written as mm/day (x86400 from kg/m2/s; 1 kg/m2 == 1 mm depth).
    // These streams use precision(4)+ios::fixed, so raw SI values ~1e-6 would all round to
    // 0.0000 — hence the unit scaling, and the _mmd suffix so the legend is unambiguous.
    // Q_precip/Q_rad are scaled to mW/m3 for the same reason. Zero unless ATJUP_PRECIP is set. This meridional cut is
    // the most useful one for precipitation: it shows the three condensation decks stacked in
    // altitude (deep H2O rain, NH4SH crystals mid-level, NH3 snow aloft) against latitude.
    // Units as above: mm/day for the fluxes, mW/m3 for the latent heating.
    dump_zonal("P_rain", P_rain, 86400.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("P_snow", P_snow, 86400.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("P_graupel", P_graupel, 86400.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("P_nh3_rain", P_nh3_rain, 86400.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("P_nh3_snow", P_nh3_snow, 86400.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("P_nh3_graupel", P_nh3_graupel, 86400.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("P_ch4_rain", P_ch4_rain, 86400.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("P_ch4_snow", P_ch4_snow, 86400.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("P_ch4_graupel", P_ch4_graupel, 86400.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("P_nh4sh", P_nh4sh, 86400.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("Q_precip", Q_precip, 1.0e3, k_zonal, Jupiter_vtk_zonal_File);

    // Turbulence (TurbulenceJup: k-epsilon / k-omega / k-omega SST). Zero unless ATJUP_TURB is
    // set. These streams use precision(4)+ios::fixed, so anything below 5e-5 rounds to 0.0000:
    // nue* is ~7e-5 and would vanish, and k* ~7e-3 would keep only two digits. The two worst-hit
    // fields are therefore written in PHYSICAL units, with the unit in the name:
    //   tke_m2s2  = k*   * u_0^2            [m2/s2]  (~67 for k* = 6.7e-3)
    //   nue_t_m2s = nue* * u_0 * L_atm[m]   [m2/s]   (~994, i.e. at ATOM's 1000 cap)
    // dis stays DIMENSIONLESS: its conversion is model-dependent (eps* uses u_0^3/L_atm, omega*
    // uses u_0/L_atm) and the writer cannot know which model ran; its raw range (3.6e-3 for
    // k-eps, 19 for k-omega) resolves fine. prod and the two source terms are left dimensionless
    // for the same reason — their magnitudes (1e-2 .. 1e2) print without loss.
    dump_zonal("tke", tke, u_0 * u_0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("dis", dis, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("nue_t", nue, u_0 * L_atm * 1.0e3, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("prod", prod, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("tke_source", tke_source, 1.0, k_zonal, Jupiter_vtk_zonal_File);
    dump_zonal("dis_source", dis_source, 1.0, k_zonal, Jupiter_vtk_zonal_File);

    Jupiter_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Jupiter_vtk_zonal_File << u.x[i][j][k_zonal] << " " << v.x[i][j][k_zonal] << " " << z << endl;
        }
    }
    ParaViewWriter<cJupiterModel>(*this).close_slice(Jupiter_vtk_zonal_File, "zonal", k_zonal, n);
    return;
}
/*
 * 
*/
void cJupiterModel::paraview_vtk_longal(int n, int j_longal){
    using namespace ParaViewIO;
    double r_mix_plus = r_mix * 1e6;
    ofstream Jupiter_vtk_longal_File = ParaViewWriter<cJupiterModel>(*this)
        .open_slice("longal", "Longitudinal", j_longal, n, km, im, 0.025, true);
    const double y = 0.0;   // out-of-plane component; longal advances z, so y stayed 0
    dump_longal("Seamount", SeaMount, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("u-Component", u, u_0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("v-Component", v, u_0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("w-Component", w, u_0, j_longal, Jupiter_vtk_longal_File);
    Jupiter_vtk_longal_File <<  "SCALARS Temperature float " << 1 << endl;
    Jupiter_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;

    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Jupiter_vtk_longal_File << t.x[i][j_longal][k] * t_ref - 273.15 << endl;
            aux.x[i][j_longal][k] = get_layer_height(i);
        }
    }

    dump_longal("thermalmassflux", thermalmassflux, 1e-3, j_longal, Jupiter_vtk_longal_File);

    dump_longal("height", aux, 1.0, j_longal, Jupiter_vtk_longal_File);

    dump_longal("CH4", ch4, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("CH4Cloud", ch4_cloud, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("CH4Ice", ch4_ice, 1.0, j_longal, Jupiter_vtk_longal_File);

    dump_longal("H2O", h2o, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("H2OCloud", h2o_cloud, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("H2OIce", h2o_ice, 1.0, j_longal, Jupiter_vtk_longal_File);

    dump_longal("H2S", h2s, 1.0, j_longal, Jupiter_vtk_longal_File);
/*
    dump_longal("w_h2s", w_h2s, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("j_h2s", j_h2s, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("jT_h2s", jT_h2s, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("massflux_h2s", massflux_h2s, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("difflux_h2s", difflux_h2s, 1.0, j_longal, Jupiter_vtk_longal_File);
*/
    dump_longal("NH3", nh3, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("NH3Cloud", nh3_cloud, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("NH3Ice", nh3_ice, 1.0, j_longal, Jupiter_vtk_longal_File);
/*
    dump_longal("w_nh3", w_nh3, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("j_nh3", j_nh3, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("jT_nh3", jT_nh3, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("massflux_nh3", massflux_nh3, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("difflux_nh3", difflux_nh3, 1.0, j_longal, Jupiter_vtk_longal_File);
*/
    dump_longal("NH4SH", nh4sh, r_mix_plus, j_longal, Jupiter_vtk_longal_File);
/*
    dump_longal("w_nh4sh", w_nh4sh, r_mix_plus, j_longal, Jupiter_vtk_longal_File);
    dump_longal("massflux_nh4sh", massflux_nh4sh, r_mix_plus, j_longal, Jupiter_vtk_longal_File);
    dump_longal("j_nh4sh", j_nh4sh, r_mix_plus, j_longal, Jupiter_vtk_longal_File);
    dump_longal("jT_nh4sh", jT_nh4sh, r_mix_plus, j_longal, Jupiter_vtk_longal_File);
    dump_longal("difflux_nh4sh", difflux_nh4sh, r_mix_plus, j_longal, Jupiter_vtk_longal_File);
*/

    dump_longal("PressureDyn", p_dyn, p_dyn_to_bar() * 1.0e3, j_longal, Jupiter_vtk_longal_File);
    dump_longal("PressureStat", p_stat, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("rho_mix", rho_mix, 1.0, j_longal, Jupiter_vtk_longal_File);

    dump_longal("CoriolisForce", CoriolisForce, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("CentrifugalForce", CentrifugalForce, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("BuoyancyForce", BuoyancyForce, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("PresGradForce", PresGradForce, 1.0, j_longal, Jupiter_vtk_longal_File);

    dump_longal("Q_Latent", Q_Latent, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("Q_Sensible", Q_Sensible, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("Q_rad_mW_m3", Q_rad, 1.0e3, j_longal, Jupiter_vtk_longal_File);
    dump_longal("Radiation", radiation, 1.0, j_longal, Jupiter_vtk_longal_File);

    // Precipitation fluxes, written as mm/day (x86400 from kg/m2/s; 1 kg/m2 == 1 mm depth).
    // These streams use precision(4)+ios::fixed, so raw SI values ~1e-6 would all round to
    // 0.0000 — hence the unit scaling, and the _mmd suffix so the legend is unambiguous.
    // Q_precip/Q_rad are scaled to mW/m3 for the same reason. Zero unless ATJUP_PRECIP is set.
    dump_longal("P_rain", P_rain, 86400.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("P_snow", P_snow, 86400.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("P_graupel", P_graupel, 86400.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("P_nh3_rain", P_nh3_rain, 86400.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("P_nh3_snow", P_nh3_snow, 86400.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("P_nh3_graupel", P_nh3_graupel, 86400.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("P_ch4_rain", P_ch4_rain, 86400.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("P_ch4_snow", P_ch4_snow, 86400.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("P_ch4_graupel", P_ch4_graupel, 86400.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("P_nh4sh", P_nh4sh, 86400.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("Q_precip", Q_precip, 1.0e3, j_longal, Jupiter_vtk_longal_File);

    // Turbulence (TurbulenceJup: k-epsilon / k-omega / k-omega SST). Zero unless ATJUP_TURB is
    // set. These streams use precision(4)+ios::fixed, so anything below 5e-5 rounds to 0.0000:
    // nue* is ~7e-5 and would vanish, and k* ~7e-3 would keep only two digits. The two worst-hit
    // fields are therefore written in PHYSICAL units, with the unit in the name:
    //   tke_m2s2  = k*   * u_0^2            [m2/s2]  (~67 for k* = 6.7e-3)
    //   nue_t_m2s = nue* * u_0 * L_atm[m]   [m2/s]   (~994, i.e. at ATOM's 1000 cap)
    // dis stays DIMENSIONLESS: its conversion is model-dependent (eps* uses u_0^3/L_atm, omega*
    // uses u_0/L_atm) and the writer cannot know which model ran; its raw range (3.6e-3 for
    // k-eps, 19 for k-omega) resolves fine. prod and the two source terms are left dimensionless
    // for the same reason — their magnitudes (1e-2 .. 1e2) print without loss.
    dump_longal("tke", tke, u_0 * u_0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("dis", dis, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("nue_t", nue, u_0 * L_atm * 1.0e3, j_longal, Jupiter_vtk_longal_File);
    dump_longal("prod", prod, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("tke_source", tke_source, 1.0, j_longal, Jupiter_vtk_longal_File);
    dump_longal("dis_source", dis_source, 1.0, j_longal, Jupiter_vtk_longal_File);

    Jupiter_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Jupiter_vtk_longal_File << u.x[i][j_longal][k] << " " 
                << y << " " << w.x[i][j_longal][k] << endl;
        }
    }
    ParaViewWriter<cJupiterModel>(*this).close_slice(Jupiter_vtk_longal_File, "longal", j_longal, n);
    return;
}
/*
 * 
*/
void cJupiterModel::paraview_sphere_vts(int n){
    using namespace ParaViewIO;
    double x, y, z, sinthe, sinphi, costhe, cosphi;
    double r_mix_plus = r_mix * 1e6;
    string Jupiter_sphere_vts_File_Name = output_path + "/Jupiter_sphere_" 
        + std::to_string(n) + ".vts";
    ofstream Jupiter_sphere_vts_File;
    Jupiter_sphere_vts_File.precision(4);
    Jupiter_sphere_vts_File.setf(ios::fixed);
    Jupiter_sphere_vts_File.open(Jupiter_sphere_vts_File_Name);
    if (!Jupiter_sphere_vts_File.is_open()){
        cerr << "ERROR: could not open paraview_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Jupiter_sphere_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Jupiter_sphere_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Jupiter_sphere_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Jupiter_sphere_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Jupiter_sphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Temperature PressureDyn PressureStat  H2S H2SCloud H2SIce CH4 CH4Cloud CH4Ice NH3 NH3Cloud NH3Ice H2O H2OCloud H2OIce \">\n"  << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        sinphi = sin( phi.z[k]);
        cosphi = cos( phi.z[k]);
        for(int j = 0; j < jm; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            // Fourth and last site of the hemispheric cosine flip; see cJupiterModel.h::
            // costhe_abs. Dead code today — the only caller of paraview_sphere_vts is
            // commented out in FileIO_Jup.cpp — but left inconsistent it would reappear
            // as a silent southern-hemisphere error the day the sphere writer is used.
            if(costhe_abs() && j > 90) costhe = - costhe;
            for(int i = 0; i < im; i++){
                // LOCALS, not aux_u/aux_v/aux_w. These are the spherical-to-Cartesian velocity
                // components and they are wanted only for the line printed immediately below —
                // but aux_* are not scratch. RHS_Jup_Turb fills aux_u with the intermediate
                // velocity rhs_u + dpdr and PressureSolver differentiates it, so writing them
                // here left the projection reading rendering coordinates.
                //
                // MEASURED ON ATSAT, where the identical loop was the whole difference: with the
                // call enabled, three iterations gave sat_restart_3.bin = a585481f... against
                // b64cb36f... without it, and a control with the panorama writing every
                // iteration and only the sphere disabled reproduced b64cb36f... exactly. This is
                // why the call site was commented out in both models rather than the writer
                // fixed — an output routine that changes the answer is worse than none.
                //
                // ATJUP's call site stays commented out for a SECOND, unrelated reason that this
                // does not address: ATJUP_METRIC_RADIUS is on by default, so rad.z runs 1..500
                // and the Cartesian construction below renders a shell of relative thickness
                // 1/500. See the note at cJupiterModel.cpp:335. Fixing that is a separate job.
                const double cart_u = sinthe * cosphi * u.x[i][j][k] + costhe * cosphi * v.x[i][j][k] - sinphi * w.x[i][j][k];
                const double cart_v = sinthe * sinphi * u.x[i][j][k] + sinphi * costhe * v.x[i][j][k] + cosphi * w.x[i][j][k];
                const double cart_w = costhe * u.x[i][j][k] - sinthe * v.x[i][j][k];
                Jupiter_sphere_vts_File << cart_u << " " << cart_v << " " << cart_w  << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << t.x[i][j][k] * t_ref - 273.15 << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureDyn\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << p_dyn.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureStat\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1e-3 * p_stat.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2O\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1.0 * h2o.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2S\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1.0 * h2s.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1.0 * nh3.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH4SH\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << r_mix_plus * nh4sh.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2OCloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1.0 * h2o_cloud.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }


    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3Cloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1.0 * nh3_cloud.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }


    Jupiter_sphere_vts_File <<  "\n"  << endl;

    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2OIce\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1.0 * h2o_ice.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }


    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3Ice\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1.0 * nh3_ice.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;


    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"CH4\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1.0 * ch4.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"CH4Cloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1.0 * ch4_cloud.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"CH4Ice\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << 1.0 * ch4_ice.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"u-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << u_0 * aux_u.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"v-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << u_0 * aux_v.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"w-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << u_0 * aux_w.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Seamount\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Jupiter_sphere_vts_File << SeaMount.x[i][j][k] << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "\n"  << endl;
    Jupiter_sphere_vts_File <<  "    </DataArray>\n"  << endl;
    Jupiter_sphere_vts_File <<  "   </PointData>\n" << endl;
    Jupiter_sphere_vts_File <<  "   <Points>\n"  << endl;
    Jupiter_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                x = rad.z[i] * sin( the.z[j])* cos(phi.z[k]);
                y = rad.z[i] * sin( the.z[j])* sin(phi.z[k]);
                z = rad.z[i] * cos( the.z[j]);
                Jupiter_sphere_vts_File << x << " " << y << " " << z  << endl;
            }
            Jupiter_sphere_vts_File <<  "\n"  << endl;
        }
        Jupiter_sphere_vts_File <<  "\n"  << endl;
    }
    Jupiter_sphere_vts_File <<  "    </DataArray>\n"  << endl;
    Jupiter_sphere_vts_File <<  "   </Points>\n"  << endl;
    Jupiter_sphere_vts_File <<  "  </Piece>\n"  << endl;
    Jupiter_sphere_vts_File <<  " </StructuredGrid>\n"  << endl;
    Jupiter_sphere_vts_File <<  "</VTKFile>\n"  << endl;
    Jupiter_sphere_vts_File.close();
    cout << "   File:  " << "Jupiter_sphere_" 
        << n << ".vts" << "  has been written to Directory:  " 
        << output_path << endl;
}
/*
 * 
*/
void cJupiterModel::JupiterPlotData(){
    // The whole body moved to the SHARED ParaViewWriter.h — it was 100 % identical between
    // the two models apart from the planet name in the file name.
    ParaViewWriter<cJupiterModel>(*this).plot_data();
}


