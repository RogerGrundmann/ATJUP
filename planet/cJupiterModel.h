#ifndef _CJUPITERMODEL_H
#define _CJUPITERMODEL_H

#include <fenv.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>    
#include <cmath>
#include <map>
#include <set>
#include <limits>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>

#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"
#include "tinyxml2.h"
#include "PythonStream.h"
#include "Utils.h"
#include "Config.h"

// Forward declarations — full definitions included at the bottom of this file
// after class cJupiterModel is complete, so inline bodies can access its members.

class ChemistryJup;
class PressureSolverJup;
class SaturationAdjustmentJup;
class BC_Jup;
class VelocityInitializerJup;


#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

namespace{
    std::function<double(double)> default_lambda=[](double i)->double{return i;};
}

class cJupiterModel{

    friend class ChemistryJup;
    friend class PressureSolverJup;
    friend class SaturationAdjustmentJup;
    friend class BC_Jup;
    friend class VelocityInitializerJup;
    friend class RadiationJup;
    friend class PrecipitationJup;
    friend class TurbulenceJup;
    friend class ConvectiveAdjustmentJup;
    friend class ThermalWindJup;

public:

    const char *filename;

    cJupiterModel();
    ~cJupiterModel();

    cJupiterModel(const cJupiterModel&) = delete;
    cJupiterModel& operator=(const cJupiterModel&) = delete;

    static cJupiterModel* get_model(){
        if(!m_model){
            m_model = new cJupiterModel();
        }
        return m_model;
    }

    void LoadConfig(const char *filename);
    void Run();

    #include "JupiterParams.h.inc"

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr;
    double dt = 0.0;
    static const double the0, phi0, r0;

    int n, n_print, panorama_cnt, iter_n;

    int i_res, j_res, k_res;

    std::vector<double> tropopause_layers; // keep the tropopause layer index
    std::vector<std::vector<int> > i_topography;
    std::vector<double> u_trans;
    std::vector<double> v_trans;
    std::vector<double> w_trans;

    double maxValue, minValue;
    /*
     * Given a latitude, return the layer index of tropopause
    */
    int get_tropopause_layer(int j){
        assert(j>=0);
        assert(j<jm);
        //refer to  BC_Thermo::TropopauseLocation and BC_Thermo::GetTropopauseHightAdd
        //tropopause height is proportional to the mean tropospheric temperature.
        //higher near the equator - warm troposphere
        //lower at the poles - cold troposphere

//        tropopause_layers[j] = 30;
//        tropopause_layers[j] = 35;
//        init_tropopause_layers();
        tropopause_layers[j] = im_tropopause[j];
        return tropopause_layers[j];
    }
    /*
     *
    */
    int get_surface_layer(int j, int k){
        return i_topography[j][k];
    }
    /*
     * This function must be called after init_layer_heights()
     * Given a layer index i, return the height of this layer
    */
    float get_layer_height(int i){
        if(0>i || i>im-1){
            return -1;
        }
        return m_layer_heights[i];
    }
    std::vector<float> get_layer_heights(){
        return m_layer_heights;
    }
    /*
     * Thickness of layer i in METRES, i.e. get_layer_height(i+1) - get_layer_height(i)
     * converted from km. NOTE: m_layer_heights is built from L_atm, which is specified in km
     * (L_atm = 140.0), so get_layer_height() returns KILOMETRES. Most callers only use it as a
     * relative height and are unaffected, but anything forming a physical per-metre quantity
     * (e.g. a volumetric heating rate W/m3 = flux difference / dz) must use this accessor —
     * using the raw km difference makes such rates 1000x too large.
    */
    double layer_thickness_m(int i){
        if(i < 0 || i > im - 2) return -1.0;
        return (double)(m_layer_heights[i + 1] - m_layer_heights[i]) * 1.0e3;
    }
    /*
    * Given a altitude, return the layer index
    */
    int get_layer_index(float height){
        std::size_t i = 0;
        for(; i<m_layer_heights.size(); i++){
            if(height<m_layer_heights[i])
                return i-1;
        }
        return i;
    }



private:

    static cJupiterModel* m_model;

    ChemistryJup* m_chem = nullptr;
    ChemistryJup& getChemistry();

    PressureSolverJup* m_pressure = nullptr;
    PressureSolverJup& getPressureSolver();

    PythonStream ps;
    std::streambuf *backup;

    const double c43 = 4.0/3.0, c13 = 1.0/3.0, c32 = 3.0/2.0, c42 = 4.0/2.0, c12 = 1.0/2.0;
    static const int im = 41, jm = 181, km = 361;

    int i_max = 40;                                                     // corresponds to about 125 km above 10e6 Pa pressure level, maximum height of the tropopause at equator
    int i_beg = 30;                                                     // corresponds to about 100 km above 10e6 Pa pressure level
//    int i_beg = 40;                                                     // corresponds to about 100 km above 10e6 Pa pressure level

    double mue_mix, k_mix, cp_mix, rg_mix, r_mix, R_mix, c_mix, x_mix;

// at 230K NH3 and H2S condense via a heterogenious reaction: NH3 + H2S -> NH4SH ( Planetary Scienses, de Pater, Lissauer)


// temperatures at triple point and ice formation
    double t_0_h2o = 273.15;                                            // in K == 0°C, triple point
    double t_00_h2o = 210.15;                                           // in K == -67°C (Planetary Siences)
//    double t_00_h2o = 241.15;                                           // in K == -32°C (COSMO)
    double t_000 = 235.15;                                              // in K == -20°C (precipitation module)

    double t_0_h2s = 187.66;                                            // in K == -53.15°C, triple point
    double t_00_h2s = 220.0;                                            // in K == -85.5°C, h2s-ice cloud formation (Planetary Sciences, p. 96)

    double t_0_nh3 = 195.5;                                             // in K == -77.65°C, triple point, gas and liquid pase
    double t_00_nh3 = 140.0;                                            // in K == -133.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)
//    double t_00_nh3 = 220.0;                                            // WRONG: 220 K > t_0_nh3=195.5 K inverts Tao mixed-phase algorithm
//    double t_00_nh3 = 120.0;                                            // in K == -153.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)

    double t_0_nh4sh = 230.0;                                           // in K == -43.15°C, nh4sh formation onset (Planetary Sciences, p. 96)
    double t_00_nh4sh = 200.0;                                          // in K == -73.15°C, nh4sh formation end (Planetary Sciences, p. 96)


// pressures take from  Planetary Sciences, p. 96
    double p_0_h2o = 2.0;                                               // in bar
    double p_00_h2o = 5.1;                                              // in bar

    double p_0_h2s = 1.3;                                               // in bar
    double p_00_h2s = 2.2;                                              // in bar

    double p_0_nh3 = 0.32;                                              // in bar
    double p_00_nh3 = 2.2;                                              // in bar

    double p_0_nh4sh = 1.3;                                             // in bar
    double p_00_nh4sh = 2.2;                                            // in bar


// constants for Clausius-Clapeyron law
    double coeff_h2_A = -2000.0;                                        // invented
    double coeff_h2_B = 8.0;                                            // invented
//    double coeff_h2_A = -125.0;                                         // AI solution
//    double coeff_h2_B = 4.0;                                            // AI solution

    double coeff_he_A = -3000.0;                                        // invented
    double coeff_he_B = 10.0;                                           // invented


    double coeff_h2o_A = -4961.04;                                      // from triple and critical point values for water
    double coeff_h2o_B = 13.0662;                                       // from triple and critical point values for water

    double coeff_h2o_A_i = -4961.04;                                    // invented
    double coeff_h2o_B_i = 13.0662;                                     // invented


    double coeff_h2s_A = -2251.66;                                      // from triple and critical point values for hydrogen sulfide
    double coeff_h2s_B = 10.5253;                                       // from triple and critical point values for hydrogen sulfide

    double coeff_h2s_A_i = -2251.66;                                    // invented
    double coeff_h2s_B_i = 10.5253;                                     // invented


    double coeff_nh3_A = -2836.56;                                      // from triple and critical point values for ammonia
    double coeff_nh3_B = 11.7271;                                       // from triple and critical point values for ammonia

    double coeff_nh3_A_i = -2836.56;                                    // invented
    double coeff_nh3_B_i = 11.7271;                                     // invented


//    double coeff_nh4sh_A = -2836.56;                                    // invented
//    double coeff_nh4sh_B = 11.7271;                                     // invented

    double coeff_nh4sh_A = -10834.0;                                    // AI
    double coeff_nh4sh_B = 34.12;                                       // AI


// molecular weights
    double m_h2 = 2.016;                                                // molecular weight of hydrogen in kg/Kmol (molar mass)
    double m_he = 4.02602;                                              // molecular weight of helium in kg/Kmol
    double m_nh3 = 17.03052;                                            // molecular weight of ammonia in kg/Kmol
    double m_nh4sh = 51.1114;                                           // molecular weight of ammonium hydrosulfide in kg/Kmol
    double m_h2s = 34.08088;                                            // molecular weight of hydrogen sulfide in kg/Kmol
    double m_h2o = 18.01588;                                            // molecular weight of water in kg/Kmol

// mass density (concentration) of species in g/cm³ == 10e6 g/m³ == 10e3 kg/m³                            Planetary sciences p. 90 2010
/*
    double rg_h2 = 0.408;                                               // mass density of vapour
    double rg_he = 0.1752;                                              // mass density of vapour
    double rg_nh3 = 0.7623;                                             // mass density of vapour
    double rg_nh4sh = 1170.0;                                           // mass density of vapour
    double rg_h2s = 1.5357;                                             // mass density of vapour
    double rg_h2o = 0.005;                                              // mass density of vapour
*/
    double rg_h2 = 0.408;                                               // mass density of vapour
    double rg_he = 0.1752;                                              // mass density of vapour
    double rg_nh3 = 817.0;                                              // AI condensed phase
    double rg_nh4sh = 1170.0;                                           // AI condensed phase
    double rg_h2s = 949.0;                                              // AI condensed phase
    double rg_h2o = 1000.0;                                             // AI condensed phase

// NH4SH crystal radius for Stokes settling
    double r_p_nh4sh = 1.0e-5;  // NH4SH crystal radius in m

// vapour mass densities of gases                Planetary sciences p. 90 2010
/*
    double r_h2 = 0.864;                                                // density of hydrogen vapour in kg/m³
    double r_he = 0.136;                                                // density of helium vapour  in kg/m³
    double r_nh3 = 0.04;                                                // density of ammonia vapour in kg/m³
    double r_h2s = 0.055;                                               // density of hydrogen sulfide vapour in kg/m³
    double r_h2o = 0.09;                                                // density of water vapour in kg/m³
    double r_nh4sh = 0.007;                                             // density of ammonium hydrosulfide vapour in kg/m³  assumption
    double r_nh3_add = 0.09;                                            // density of ammonia vapour in kg/m³
*/
    double r_h2 = 0.864;                                                // density of hydrogen vapour in kg/m³
    double r_he = 0.136;                                                // density of helium vapour  in kg/m³
    double r_nh3 = 0.004;                                               // density of ammonia vapour in kg/m³
    double r_h2s = 0.0004;                                              // density of hydrogen sulfide vapour in kg/m³
    double r_h2o = 0.09;                                                // density of water vapour in kg/m³
    double r_nh4sh = 0.000;                                             // density of ammonium hydrosulfide vapour in kg/m³  assumption

    // Deep, WELL-MIXED ammonia below the condensation level. init_vapour_cloud_ice uses it for
    // every cell warmer than t_add (= t_0_nh3), and the saturation profile capped at r_nh3
    // everywhere above; so this value is the deep plateau and r_nh3 is the ceiling aloft.
    //
    // It was 0.0001, i.e. FORTY TIMES BELOW r_nh3 = 0.004, which inverted the profile: NH3 sat
    // at 8.3e-5 from i=0 to i=12 and then rose 40x to the 3.5e-3 cap at i=22 (77 km), a vapour
    // MAXIMUM above the cloud deck. On Jupiter NH3 is well mixed at depth and decreases upward
    // through the NH4SH and NH3 decks, so the profile ran the wrong way round. The consequence
    // was measurable: NH4SH forms exactly at the crossing, i=12..13, where the 40x gradient puts
    // the whole ammonia reservoir in one pair of layers, and the field peaked at 0.33 g/m3 there
    // with NH3 locally drawn down from 8.5e-5 to 1.3e-5.
    //
    // The reference block above this one (Planetary Sciences p. 90) has the relation the right
    // way round, r_nh3_add = 0.09 against r_nh3 = 0.04, a factor 2.25. That ratio is carried over
    // to the abundance scale of the active block rather than importing the absolute reference
    // values, which would change far more than the profile shape. With r_nh3_add >= r_nh3 the
    // initial profile is monotonically non-increasing upward by construction: constant
    // r_nh3_add while t > t_add, then min(saturation, r_nh3) <= r_nh3 <= r_nh3_add.
    //
    // Override at runtime with ATJUP_R_NH3_ADD [kg/m3] to A/B the deep abundance without a
    // rebuild; 0 or unset keeps the value below.
    double r_nh3_add = [](){
        const char* e = getenv("ATJUP_R_NH3_ADD");
        const double v = e ? atof(e) : 0.0;
        return (v > 0.0) ? v : 0.009;
    }();                                                                // density of ammonia vapour in kg/m³

// vapour mass densities of clouds and ices               Planetary sciences p. 90 2010
    double r_nh3_cloud = 0.004;                                         // density of ammonia cloud in kg/m³
    double r_h2o_cloud = 0.009;                                         // density of water cloud in kg/m³
    double r_nh3_ice = 0.0003;                                          // density of ammonia ice in kg/m³
    double r_h2o_ice = 0.002;                                           // density of water ice in kg/m³

// vapour molar densities of gases
    double c_h2 = r_h2/m_h2;                                            // density of hydrogen vapour in kmol/m³
    double c_he = r_he/m_he;                                            // density of helium vapour  in kmol/m³
    double c_nh3 = r_nh3/m_nh3;                                         // density of ammonia vapour in kmol/m³
    double c_h2s = r_h2s/m_h2s;                                         // density of hydrogen sulfide vapour in kmol/m³
    double c_h2o = r_h2o/m_h2o;                                         // density of water vapour in kmol/m³
    double c_nh4sh = r_nh4sh/m_nh4sh;                                   // density of ammonium hydrosulfide vapour in kmol/m³  assumption
/*
 // ratio of vapour molecular weight to mean molecular weight
    double X_h2 = 0.864;
    double X_he = 0.136;
    double X_h2o = 1.7e-3;
    double X_nh3 = 2.0e-4;
    double X_h2s = 7.7e-5;
    double X_nh4sh = 3.6e-5;
*/
 // mass fraction
    double X_h2 = 0.864;
    double X_he = 0.136;
    double X_h2o = 1.7e-3;                                              // AI
    double X_nh3 = 2.0e-4;                                              // AI
    double X_h2s = 7.7e-5;                                              // AI
    double X_nh4sh = 0.0;                                               // AI

// gas constants
    double R_h2 = 4124.2;                                               // gas constant of hydrogen in J/(kg*K)
    double R_he = 2077.1;                                               // gas constant of helium in J/(kg*K)
    double R_nh3 = 488.21;                                              // gas constant of ammoinia in J/(kg*K)
    double R_nh4sh = 162.67;                                            // gas constant of ammoinium hydrosulfide in J/(kg*K)   by AI
    double R_h2s = 243.96;                                              // gas constant of hydrogen sulfide inJ/(kg*K) 
    double R_h2o = 461.52;                                              // gas constant of water inJ/(kg*K)
/*
// dynamic viscosities
    double mue_h2 = 0.84e-5;                                            // dynamic viscosity of hydrogen in Ns/m²
    double mue_he = 1.87e-5;                                            // dynamic viscosity of helium in Ns/m²
    double mue_nh3 = 0.92e-5;                                           // dynamic viscosity of ammonia in Ns/m²
    double mue_nh4sh = 0.99e-5;                                         // dynamic viscosity of ammonium sulfide in Ns/m²
    double mue_h2s = 1.3e-5;                                            // dynamic viscosity of hydrogen sulfid in Ns/m²
    double mue_h2o = 1.308e-3;                                          // dynamic viscosity of water in Ns/m²
*/
// dynamic viscosities
    double mue_h2 = 0.84e-5;                                            // dynamic viscosity of hydrogen in Ns/m²
    double mue_he = 1.87e-5;                                            // dynamic viscosity of helium in Ns/m²
    double mue_nh3 = 0.92e-5;                                           // dynamic viscosity of ammonia in Ns/m²
    double mue_nh4sh = 0.0;                                             // dynamic viscosity of ammonium sulfide in Ns/m²
    double mue_h2s = 1.3e-5;                                            // dynamic viscosity of hydrogen sulfid in Ns/m²
    double mue_h2o = 0.9e-5;                                            // dynamic viscosity of water in Ns/m²   by AI

// thermal conductivities
    double k_h2 = 0.1317;                                               // thermal conductivity of hydrogen in W/(m*K)
    double k_he = 0.1193;                                               // thermal conductivity of helium in W/(m*K)
    double k_nh3 = 0.02102;                                             // thermal conductivity of ammonia in W/(m*K)
    double k_nh4sh = 0.0;                                               // thermal conductivity of ammonium hydrosulfide in W/(m*K)
    double k_h2s = 0.013;                                               // thermal conductivity of hydrogen sulfide in W/(m*K)
    double k_h2o = 0.0187;                                              // thermal conductivity of water in W/(m*K)

// specific heat capacities
    double cp_h2 = 14.32e3;                                             // specific heat capacity of hydrogen in J/(kg*K)
    double cp_he = 5.19e3;                                              // specific heat capacity of helium in J/(kg*K)
    double cp_nh3 = 2.19e3;                                             // specific heat capacity of ammonia in J/(kg*K)
    double cp_nh4sh = 2.00e3;                                           // specific heat capacity of ammonium hydrosulfide in J/(kg*K)
    double cp_h2s = 2.24e3;                                             // specific heat capacity of hydrogen sulfid in J/(kg*K)
    double cp_h2o = 1.93e3;                                             // specific heat capacity of water in J/(kg*K)
/*
// ratios of gas constants of dry gas to vapour or vapour molecular weight to mean atmospheric molecular weight or m/m_mix
    double ep_h2 = 0.8572;                                              // ratio of the gas constants of dry air to h2 non-dimensional or m/m_mix
    double ep_he = 1.7152;                                              // ratio of the gas constants of dry he to h2 non-dimensional
    double ep_h2o = 8.1253;                                             // ratio of the gas constants of dry air to h2 non-dimensional
    double ep_h2s = 14.5192;                                            // ratio of the gas constants of dry hydrogen sulfide to h2 non-dimensional
    double ep_nh3 = 7.6752;                                             // ratio of the gas constants of dry ammonia to h2 non-dimensional
    double ep_nh4sh = 21.7745;                                          // ratio of the gas constants of dry ammonia hydrosufide to h2 non-dimensional       invented
*/
// ratios of gas constants of dry gas to vapour or vapour molecular weight to mean atmospheric molecular weight or m/m_mix
    double ep_h2 = 0.908;                                               // ratio of the gas constants of dry air to h2 non-dimensional or m/m_mix
    double ep_he = 1.803;                                               // ratio of the gas constants of dry he to h2 non-dimensional
    double ep_h2o = 8.115;                                              // ratio of the gas constants of dry air to h2 non-dimensional
    double ep_h2s = 15.35;                                              // ratio of the gas constants of dry hydrogen sulfide to h2 non-dimensional
    double ep_nh3 = 7.671;                                              // ratio of the gas constants of dry ammonia to h2 non-dimensional
    double ep_nh4sh = 23.02;                                            // ratio of the gas constants of dry ammonia hydrosufide to h2 non-dimensional       invented

// latent heat of evaporation
    double lv_h2o = 2.5009e6;                                           // latent heat of h2o evaporation at 0°C in J/kg
    double lv_h2s = 3.5340e6;                                           // latent heat of h2s evaporation at -73°C in J/kg
    double lv_nh3 = 1.3720e6;                                           // latent heat of nh3 evaporation at -33.33 in J/kg

// latent heat of sublimation
    double ls_h2o = 2.8339e6;                                           // latent heat of h2o sublimation at 0°C in J/kg
    double ls_h2s = 7.4500e5;                                           // latent heat of h2s sublimation at -98°C in J/kg
    double ls_nh3 = 1.8320e6;                                           // latent heat of nh3 sublimation at -93.15 in J/kg
    double ls_nh4sh = 1.7e6;                                            // latent heat of nh3 sublimation at -93.15 in J/kg

// Schmidt number
    double sc_h2 = 0.20;                                                // Schmidt numbert of h2o, Sc = nue/D
    double sc_he = 0.22;                                                // Schmidt numbert of h2o, Sc = nue/D
    double sc_h2o = 0.61;                                               // Schmidt numbert of h2o, Sc = nue/D
    double sc_h2s = 0.94;                                               // Schmidt number of h2s, Sc = nue/D 
    double sc_nh3 = 0.61;                                               // Schmidt number of nh3, Sc = nue/D 
    double sc_nh4sh = 0.7;                                              // Schmidt number of nh4sh, Sc = nue/D 

// Prantl numbers
    double Pr = 0.72;                                                   // Prandtl-number 
/*
// diffusion coefficients
    double D_nh3 = 1.5e-9;                                              // ordinary diffusion coefficient of ammonia in m*m/s 
    double D_h2s = 1.36e-9;                                             // ordinary diffusion coefficient of hydrogen sulfid in m*m/s 
    double D_nh4sh = 1.45e-9;                                           // ordinary diffusion coefficient of ammonium hydrosulfide in m*m/s
*/
// diffusion coefficients
    double D_nh3 = 1.0e-4;                                              // ordinary diffusion coefficient of ammonia in m*m/s 
    double D_h2s = 1.0e-4;                                              // ordinary diffusion coefficient of hydrogen sulfid in m*m/s 
    double D_nh4sh = 1.0e-4;                                            // ordinary diffusion coefficient of ammonium hydrosulfide in m*m/s
/*
// thermal diffusion coefficients
    double DT_nh3 = 1.54e-9;                                            // thermal diffusion coefficient of ammonia in kg/(s*m)     unklar
    double DT_h2s = 1.36e-9;                                            // thermal diffusion coefficient of hydrogen sulfid in kg/(s*m)
    double DT_nh4sh = 1.45e-9;                                          // thermal diffusion coefficient of ammonium hydrosulfide in kg/(s*m)
*/
// thermal diffusion coefficients
    double DT_nh3 = 0.0;                                                // thermal diffusion coefficient of ammonia in kg/(s*m)     unklar
    double DT_h2s = 0.0;                                                // thermal diffusion coefficient of hydrogen sulfid in kg/(s*m)
    double DT_nh4sh = 0.0;                                              // thermal diffusion coefficient of ammonium hydrosulfide in kg/(s*m)

// constants for saturation vapour pressure and latent heat from the original paper by Sanchez-Lavega, Perez-Hoyos and Huesco, p. 770
    double C_h2o = 25.096;                                              //  in bar
    double C_nh3 = 27.863;                                              //  in bar
    double C_h2s = 17.064;                                              //  in bar
    double C_nh4sh = 75.678;                                            //  in bar

    double L0_h2o = 3148.2;                                             //  in J/g
    double L0_nh3 = 2016.0;                                             //  in J/g
    double L0_h2s = 747.0;                                              //  in J/g
    double L0_nh4sh = 2915.7;                                           //  in J/g

// alf and bet are empirical constants for each phase
    double del_alf_h2o = 0.0;
    double del_bet_h2o = - 8.7e-3;

    double del_alf_nh3 = - 0.888;                                       // original paper by Sanchez-Lavega, Perez-Hoyos and Huesco, p. 770
    double del_bet_nh3 = 0.0;

// ice-phase (sublimation) Sanchez-Lavega SVP parameters — derived so that:
//   E_ice(T_triple) == E_liquid(T_triple)  (continuity at triple point)
//   L0_ice = L0_liquid * (ls / lv)         (sublimation latent heat scales the slope)
//   del_alf / del_bet kept equal to liquid  (same temperature-dependence structure)
// Result: E_ice < E_liquid for all T < T_triple, as required by Tao mixed-phase scheme.
    double C_h2o_ice    = 28.418;                                       // H2O ice SVP constant [bar]; calibrated at T_tp=273.16 K
    double L0_h2o_ice   = 3567.3;                                       // H2O L0_ice = 3148.2*(ls/lv) = 3148.2*(2833.9/2500.9) [J/g]
    double del_alf_h2o_ice = 0.0;
    double del_bet_h2o_ice = -8.7e-3;

    double C_nh3_ice    = 34.948;                                       // NH3 ice SVP constant [bar]; calibrated at T_tp=195.4 K
    double L0_nh3_ice   = 2692.5;                                       // NH3 L0_ice = 2016.0*(ls/lv) = 2016.0*(1832/1372) [J/g]
    double del_alf_nh3_ice = -1.2;
    double del_bet_nh3_ice = 0.0;

    double del_alf_h2s = 0.0;
    double del_bet_h2s = - 2.9e-3;

    double del_alf_nh4sh = - 1.760;
    double del_bet_nh4sh = 7.8e-4;

// ============================================================================
// CH4 (methane) parameters — ported from ATNEPT
// ============================================================================
    double t_0_ch4  = 90.69;                                            // in K, triple point
    // WAS 190.56 K, which is methane's CRITICAL temperature, not an ice-cloud bound — it sat 100 K
    // ABOVE t_0_ch4 and so inverted the ordering every other species has. The consequence was exact
    // and total: PrecipitationJup/Sat gate ice autoconversion on (T < t_frz && T >= t_low), which
    // for CH4 read (T < 90.69 && T >= 190.56) — an EMPTY band, so methane ice could never convert
    // to snow at any temperature. ATSAT was carrying a 45.6 g/m3 methane ice deck with zero methane
    // snow, and its clamp budget showed ch4_ice clipping 110 % of its own mass per 12 iterations
    // because the field had no sink at all.
    //
    // 67.36 K is t_0_ch4 scaled by the H2O/NH3 proportion, on Roger's instruction: H2O sits at
    // 210.15/273.15 = 0.7694 of its triple point and NH3 at 140.0/195.5 = 0.7161, mean 0.7427, and
    // 0.7427 * 90.69 = 67.36. That makes the CH4 ice band 67.36 .. 90.69 K. It is a proportion
    // carried across from two other substances, not a measured property of methane ice.
    double t_00_ch4 = 67.36;   // in K == -205.79 degC, ch4-ice cloud formation
    double p_0_ch4  = 0.1;                                              // in bar
    double p_00_ch4 = 1.1;                                              // in bar

    double coeff_ch4_A   = -1033.3;                                     // Clausius-Clapeyron, methane
    double coeff_ch4_B   = 6.3910;
    double coeff_ch4_A_i = -1033.3;                                     // ice phase, invented
    double coeff_ch4_B_i = 6.3910;

    double rho_cond_ch4 = 0.657;                                        // liquid methane density [kg/m³]
    double rg_ch4       = 0.657;                                        // mass density of ch4 vapour
    double m_ch4        = 16.042;                                       // molecular weight of methane [kg/kmol]
    double r_ch4        = 0.19;                                         // density of ch4 vapour [kg/m³]
    double r_ch4_ice    = 0.0082;                                       // density of ch4 ice [kg/m³]
    double c_ch4        = r_ch4 / m_ch4;                                // molar density of ch4 vapour [kmol/m³]
    double X_ch4        = 3.6e-5;                                       // mass fraction
    double R_ch4        = 518.28;                                       // gas constant of methane [J/(kg*K)]
    // 1.107e-5 Pa*s, not 1.107e-2. Methane's GAS viscosity is 1.107e-2 CENTIPOISE, and the
    // centipoise figure was entered as if it were Pa*s — a factor 1000. It mattered because
    // ThermalPropertiesJup forms mue_mix as a mass-weighted mean and r_ch4 = 0.19 is the third
    // largest mass fraction, so methane alone contributed 2.10e-3 of a 2.11e-3 sum: mue_mix came
    // out 1.646e-3 Ns/m2 where every component is of order 1e-5. The only consumer is the NH4SH
    // Stokes settling velocity in rhs_nh4sh, which goes as 1/mue_mix and was therefore 166x too
    // slow — 4.1e-4 m/s against the 7.7e-2 m/s that PrecipitationJup computes for the same
    // crystals from its own hand-entered viscosity.
    double mue_ch4      = 1.107e-5;                                     // dynamic viscosity of methane [Ns/m²]
    double k_ch4        = 0.0;                                          // thermal conductivity of methane [W/(m*K)]
    double cp_ch4       = 2.232e3;                                      // specific heat capacity of methane [J/(kg*K)]
    double ep_ch4       = 7.6752;                                       // R_h2/R_ch4 ratio
    double lv_ch4       = 5.11e5;                                       // latent heat of evaporation [J/kg]
    double ls_ch4       = 5.11e5;                                       // latent heat of sublimation [J/kg]
    double sc_ch4       = 0.99;                                         // Schmidt number of ch4
    double D_ch4        = 1.0e-4;                                       // ordinary diffusion coefficient [m²/s]
    double DT_ch4       = 0.0;                                          // thermal diffusion coefficient [kg/(s*m)]
    double C_ch4        = 1.627;                                        // SVP constant [bar]
    double L0_ch4       = 553.1;                                        // SVP slope constant [J/g]
    double del_alf_ch4  = 1.002;
    double del_bet_ch4  = -4.1e-3;
    // ATNEPT only specifies one CH4 SVP curve (ls_ch4 == lv_ch4); reuse the same
    // parameters for the ice-phase calibration that SaturationAdjustmentJup expects.
    double C_ch4_ice       = 1.627;
    double L0_ch4_ice      = 553.1;
    double del_alf_ch4_ice = 1.002;
    double del_bet_ch4_ice = -4.1e-3;

    std::vector<std::vector<int> > j_ellipse;
    bool has_welcome_msg_printed;
    double out_maxValue() const;
    double out_minValue() const;

    void init_layer_heights(){
        float h = L_atm/(im-1);
        for(int i=0; i<im; i++){
            m_layer_heights.push_back(i * h);
        }
        return;
    }

    // Jupiter's mean radius [km]. Deliberately a named constant rather than a config parameter:
    // nothing in the model reads a planetary radius today, which is the whole problem. It exists
    // so checkMetricConsistency() has something to compare the metric against, and so
    // ATJUP_METRIC_RADIUS has a documented target instead of a number people copy from a note.
    static constexpr double R_planet_km = 69911.0;

    // Physical length that ONE unit of rad.z represents, in km.
    //
    // Simpler here than in ATOM, and worth stating because there it is the trap: ATJUP's
    // init_layer_heights above is LINEAR (i * L_atm/(im-1)) — the exponentially stretched variant
    // sits commented out right below it — and the shell spans exactly (im-1)*dr = 1.0 in rad.z.
    // So one rad.z unit is L_atm, 140 km, full stop. In ATOM the same quantity is the shell
    // thickness while L_atm is only the stretch amplitude, and the two differ by 1/dr = 40.
    double metricShellLength_km() const {
        const double span = rad.z[im-1] - rad.z[0];
        return (span > 0.0) ? (L_atm / span) : L_atm;
    }

    // Startup check: does the radius the metric uses equal the radius the planet has?
    void checkMetricConsistency() const;

/*
    void init_layer_heights(){
        const float zeta = 3.715;
        float h = L_atm;
        for(int i=0; i<im; i++){
            m_layer_heights.push_back((exp(zeta 
            * (rad.z[i] - 1.0)) - 1) * h);                              // in m      local atmospheric shell thickness
//            std::cout << m_layer_heights.back() << std::endl;
        } 
        return;
    }
*/

    struct CellGeometry {
        double rm, rm2, exp_rm, exp_2_rm;
        double sinthe, sinthe2, costhe, cotanthe;
        double inv_rm, inv_rm2;
        double inv_rmsinthe, inv_rm2sinthe, inv_rm2sinthe2;
        double costhe_inv_rm2sinthe;
        double inv_2dr, inv_2dthe, inv_2dphi;
        double inv_dr2, inv_dthe2, inv_dphi2;
    };

    void SetDefaultConfig();
    void RHSJup(int i, int j, int k, const CellGeometry& geo);
    void RungeKuttaJup();

    // Shapiro de-checkerboarding of the velocity fields (u,v,w). Order/strength are
    // read once from the environment (see cJupiterModel.cpp):
    //   ATJUP_VEL_SHAPIRO_ORDER    2 = 1-2-1 (default, bit-identical), 4 = shear-preserving
    //   ATJUP_SHAPIRO_STRENGTH     filter strength (default 1.0)
    //   ATJUP_VEL_SHAPIRO_INLOOP   # of per-iteration passes after RK4 (default 0 = off)
    // dampVelocities() applies the configured init-time filter; the in-loop knob is
    // consulted directly in Run().
    void dampVelocities();
    void JupiterPlotData();
    void paraview_vtk_longal(int n, int j_longal);
    void paraview_vtk_radial(int n, int i_radial);
    void paraview_vtk_zonal(int n, int k_zonal);
    void paraview_panorama_vts(int n);
    void paraview_sphere_vts(int n);

    void searchMinMax_3D(string, string, 
        string, Array &, double coeff=1., 
        std::function< double(double) > lambda = default_lambda,
        bool print_heading=false);

    void searchMinMax_2D(string, string, 
        string, Array_2D &, double coeff=1.0);

    void print_welcome_msg();
    void print_final_msg();
    void printMinMax();
    void initMsg();
    void writeResults();
    void writeData();

    // Boussinesq base state for the buoyancy term: the area-weighted horizontal mean of
    // ATJUP's own buoyancy expression at each radial level, over fluid cells only. Refilled
    // once per RK4 step. See RungeKutta_Jup_Turb.cpp for why this exists.
    std::vector<double> buoy_ref_level;
    void computeBuoyancyRefLevel();
    void computeHydrostaticPressure();

    // Fills wall_nue from the SeaMount contour. Call once, after bcSeaMount().
    void computeWallViscosity();

    // Local mixture density rho = p/(R_mix*T); call with Forces(). See InitValues_Jup.cpp.
    void computeMixtureDensity();

    // ATJUP_LOCAL_RHO=1 — use the LOCAL mixture density rho_mix in place of the constant
    // reference density r_mix wherever a density genuinely belongs to the cell: the saturation
    // vapour density in SaturationAdjustmentJup and PrecipitationJup, and the Lewis groups in
    // ChemistryJup. Default 0 keeps r_mix everywhere and is bit-identical.
    //
    // Why it is a knob and not simply a fix: r_mix = 1.2844 kg/m3 is the density of the DEEPEST
    // layer, and the shell spans 0.006 to 1.09 kg/m3 — nearly three orders of magnitude. So this
    // is not a small correction, it rescales the saturation vapour density by up to 1/200 near
    // the top, and (rho/r_mix)^2 in the thermo-diffusion flux jT_*. Measure before adopting.
    // cos(theta) must change SIGN across the equator. Four places used to force it positive
    // with `if(j > 90) costhe = -costhe`, which is wrong under any velocity-sign convention:
    // the horizontal Coriolis parameter f = 2*Omega*cos(theta) would then have the same sign in
    // both hemispheres — cyclones turning the same way north and south — and cot(theta), which
    // rides on the same cosine, would be wrong in the south wherever it appears: the d/dtheta
    // term of every Laplacian (costhe_inv_rm2sinthe) and the curvature groups of the transport
    // terms. The codebase was already split on it: PressureSolverJup and ChemistryJup use the
    // unflipped cosine, so the two halves of the model disagreed south of the equator.
    //
    // the.z runs 0 (north pole) to pi (south pole), so cos(the.z[j]) is already the correct
    // colatitude cosine and the flip only destroyed it. ATJUP_COSTHE_ABS=1 restores the old
    // behaviour for A/B work.
    // Polar metric floor: sin(theta) is held at this value instead of going to zero, so that
    // 1/sin(theta) and 1/sin^2(theta) in the metric terms stay bounded. 0.55 corresponds to
    // theta = 33.4 deg, i.e. the floor is active POLEWARD OF 56.6 DEGREES LATITUDE — 16.5 % of
    // the sphere, both polar caps, where every metric term is consequently wrong.
    //
    // It is a genuine stability measure, not an oversight: the sequential k-loop of the Runge-
    // Kutta makes an asymmetric phi-Laplacian whose error is amplified by 1/sin^2, and the value
    // was raised from 0.4 to 0.55 to stop a long-run polar blow-up. Lowering it is therefore a
    // trade, not a fix, which is why it stays at 0.55 by default and is now merely VISIBLE and
    // measurable rather than hard-coded in two files that had to be kept in step by hand.
    // ATJUP_SINTHE_MIN=<x> sets it; see the measurement in the commit that introduced the knob.
    static double sinthe_min(){
        static const double v = [](){
            const char* e = getenv("ATJUP_SINTHE_MIN");
            const double x = e ? atof(e) : 0.55;
            return (x > 0.0 && x < 1.0) ? x : 0.55;
        }();
        return v;
    }

    static bool costhe_abs(){
        static const bool v = [](){ const char* e = getenv("ATJUP_COSTHE_ABS"); return e && atoi(e) != 0; }();
        return v;
    }

    static bool local_rho(){
        static const bool v = [](){ const char* e = getenv("ATJUP_LOCAL_RHO"); return e && atoi(e) != 0; }();
        return v;
    }

    // ---- What p_dyn actually is, and what it is worth in bar ----
    //
    // p_dyn is not a pressure in bar, although five places read it as one. Nothing in the model
    // ever assigns it a physical pressure: it is created solely by PressureSolverJup, which
    // relaxes  lap(p_dyn) = div(aux),  and aux_u/v/w are the momentum RHS without their own
    // pressure gradient — that is, NONDIMENSIONAL accelerations in units of u_0^2/L_atm. The
    // units of p_dyn follow from that equation and from nothing else: if lap has 1/length^2 and
    // div has 1/length, then p_dyn carries one power of length more than aux, so grad(p_dyn) is
    // an acceleration in the same units as aux — which is exactly how RHS_Jup_Turb.cpp uses it.
    // Written out, p_dyn = p_phys / (rho * u_0^2): the nondimensional kinematic pressure.
    //
    // TWO CONSEQUENCES, and they pull in opposite directions from what the earlier note in
    // RHS_Jup_Turb.cpp assumed.
    //
    // 1. The pressure gradient in the momentum equation needs NO conversion factor. The 7.79 =
    //    1e5/(r_mix*u_0^2) that ATJUP_PGRAD_SCALE offers is derived from p_dyn being in bar, and
    //    it is not. The projection is already self-consistent — grad(p_dyn) is nondimensional
    //    because the Poisson equation made it so — which is also why raising ATJUP_PGRAD_SCALE
    //    on its own measured WORSE rather than better: it was not restoring a balance, it was
    //    breaking one. The knob stays for the record, at 1.0, and 1.0 is the correct value.
    //
    // 2. Everywhere p_dyn is ADDED to p_stat, which really is in bar, it must be converted
    //    first. That is this factor:  p_bar = p_dyn * r_mix * u_0^2 / 1e5 = p_dyn * 0.128.
    //    Unconverted, p_dyn entered the density and the buoyancy 7.79x too strongly. The
    //    affected places are the buoyancy anomaly (RHS_Jup_Turb.cpp), its level mean
    //    (computeBuoyancyRefLevel), the mixture density and Forces() (InitValues_Jup.cpp) and
    //    the printed/exported field (PrintMsg_Jup.cpp, ParaView_Jup.cpp).
    //
    // ATJUP_PDYN_UNITS=0 restores the old reading (p_dyn taken as bar) for A/B work.
    static bool pdyn_units(){
        static const bool v = [](){ const char* e = getenv("ATJUP_PDYN_UNITS"); return e ? atoi(e) != 0 : true; }();
        return v;
    }
    // Multiply a p_dyn value by this to get bar. 1.0 with the knob off, i.e. the old behaviour.
    double p_dyn_to_bar() const {
        return pdyn_units() ? r_mix * u_0 * u_0 * 1.0e-5 : 1.0;
    }

    // The pressure the BUOYANCY density is built from, in bar. Hydrostatic only: p_dyn is a
    // Lagrange multiplier for the velocity constraint, not a thermodynamic variable, and putting
    // it in the equation of state closes a feedback loop that destroys the run as soon as the
    // pressure equation is solved rather than smeared. See the long note in RHS_Jup_Turb.cpp.
    // ATJUP_BUOY_PDYN=1 restores the old (p_stat + p_dyn) reading.
    static bool buoy_pdyn(){
        static const bool v = [](){ const char* e = getenv("ATJUP_BUOY_PDYN"); return e && atoi(e) != 0; }();
        return v;
    }
    double buoy_pressure(int i, int j, int k) const {
        return buoy_pdyn() ? p_stat.x[i][j][k] + p_dyn.x[i][j][k] * p_dyn_to_bar()
                           : p_stat.x[i][j][k];
    }
    // The density to use in such a place. Falls back to r_mix where rho_mix is not positive
    // (solid cells, or before the first computeMixtureDensity), so no caller can divide by zero.
    double rho_at(int i, int j, int k) const {
        if(!local_rho()) return r_mix;
        const double rho = rho_mix.x[i][j][k];
        return (rho > 0.0 && std::isfinite(rho)) ? rho : r_mix;
    }

    // Binary checkpoint / restart of the full 3D state (FileIO_Jup.cpp), the ATJUP
    // counterpart of ATOM's cAtmosphereModel::save_state / load_state.
    std::vector<Array*> restart_arrays();   // the prognostic 3D fields a checkpoint serializes
    void save_state(int iter);              // dump restart_arrays() to output_path/jup_restart_<iter>.bin
    bool load_state(int iter);              // restore them; false (and run from scratch) if absent/mismatched
    bool restart_state_is_clean();          // true when every serialized field is finite everywhere
    bool nan_watch(int iter);               // ATJUP_NANCHECK: report the first non-finite cell, by field name
    void momentum_profile(int iter);        // ATJUP_WPROFILE: per-radial-level census of u,v,w

    // Zero floor for the condensable species, with accounting — see FileIO_Jup.cpp.
    void clampNegativeSpecies();
    void reportClampBudget();
    std::vector<double> clamp_added;        // cumulative mass added by the floor, per field
    std::vector<long>   clamp_cells;        // cumulative number of clipped cells, per field

    void BC_phi();
    void BC_radius();
    void BC_theta();
    void BC_seamount();
    void BC_solidground();
    void BC_vel_surf_sur();
    void BC_scalar_surf_sur();
    void resetArrays();
    void TropopauseLocation();
    void JupiterCellStructure();
    void JupiterCellStructure_new();
    void Jupiter_PlotData();

    void init_tropopause_layers();

//    void computePressure();
    void init_temperature();
    void init_PressureStatic();
    void init_PressureDynamic();
//    void init_Density();

    void init_vapour(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, double &X,
        Array &c, Array &cloud, Array &ice);

    void init_vapour_cloud_ice(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet,
        Array &c, Array &cloud, Array &ice);

    void init_h2s(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, Array &c);

    void init_nh4sh(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, 
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, Array &c);

    void Saturation_Adjustment(std::string gas, 
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00, 
        double &ep, double &lv, double &ls, double &cp, double &r,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, double &m,
        Array &c, Array &cloud, Array &ice);

    void OneCategoryIceScheme();

    void ChemMassRateJup();
    void DiffMassFluxJup();
    void ThermalPropertiesJup();
    void Latent_Heat();
    void Forces();

    void steadyQuery();
    void restoreVar(double coeff);

    double saturation_vapour_pressure(double &T_K, double &C,
        double &L0, double &R,
        double &del_alf, double &del_bet);

    double Clausius_Clapeyron(double &T_K, double &A, double &B);  // temperature in °K
    double Humility_critical(double &x, double Hu_cr_max, double Hu_cr_mid);



    std::vector<int> im_tropopause; // keep the tropopause layer index
    std::vector<float> m_layer_heights;
    std::vector<double> cloud_loc; // lateral cloudwater distribution
    std::vector<double> r_max; // lateral r_max distribution
    std::vector<double> r_max_add; // lateral r_max distribution
    std::vector<double> t_add; // lateral r_max distribution

    Array_1D rad;
    Array_1D the;
    Array_1D phi;

    // First isothermal ("stratospheric") layer of each column, i.e. the lowest i at which the
    // tropospheric lapse rate was clamped by init_temperature; im when the column never clamps.
    // init_PressureStatic needs it because p_stat ~ T^(g/(gam*R_ref)) is the POLYTROPIC relation
    // and holds only where the lapse rate really is gam — above the clamp the pressure must
    // follow the isothermal hydrostatic law instead.
    std::vector<std::vector<int> > i_strato_base;

    // Snapshot of the lid temperature t.x[im-1][j][k] taken from the initial condition on the
    // first bcRadius() call (ATOM's t_top_init). Empty until then; used only when the
    // ATJUP_BC_T_LID_PIN knob is set. See the discussion in BC_Jup.h bcRadius().
    std::vector<std::vector<double> > t_top_init;

    Array_2D Topography; // topography
    Array_2D LatentHeat;        // areas of higher latent heat
    Array_2D Precipitation;        // areas of higher precipitation
    // All-species surface precipitation map [kg/m2/s], filled by PrecipitationJup at the base of
    // each column (i_topography, so the GRS solid is respected). Per-species totals plus their
    // sum, so a single lat-lon view shows which condensate dominates where. Zero unless
    // ATJUP_PRECIP is set. ParaView writes these as mm/day.
    Array_2D precip_srf_total;  // H2O + NH3 + CH4 + NH4SH
    Array_2D precip_srf_h2o;    // H2O rain + snow + graupel
    Array_2D precip_srf_nh3;    // NH3 rain + snow + graupel
    Array_2D precip_srf_ch4;    // CH4 rain + snow + graupel
    Array_2D precip_srf_nh4sh;  // NH4SH settling crystals
    Array_2D precipitable_water;// areas of precipitable water in the air
    Array_2D nh3_total;            // areas of higher nh3 concentration
    Array_2D nh3_cloud_total;    // areas of higher nh3_cloud concentration
    Array_2D nh3_ice_total;        // areas of higher nh3_ice concentration
    Array_2D aux_2D_v;            // auxilliar field v
    Array_2D aux_2D_w;            // auxilliar field w
    Array_2D tropopause_height; // local height of the tropopause

    Array t;                    // temperature
    Array u;                    // u-component velocity component in r-direction
    Array v;                    // v-component velocity component in theta-direction
    Array w;                    // w-component velocity component in phi-direction

    Array rho_mix;                    // density of mixture

    Array h2o;                    // water vapour
    Array h2o_cloud;                // cloud water
    Array h2o_ice;                    // cloud ice
    Array h2s;                    // water vapour
    Array nh3;                    // nh3-vapour
    Array nh3_cloud;            // nh3-cloud
    Array nh3_ice;                // nh3-ice
    Array ch4;                    // ch4-vapour
    Array ch4_cloud;            // ch4-cloud
    Array ch4_ice;                // ch4-ice
    Array nh4sh;                    // nh4sh-vapour

    Array tn;                    // temperature new
    Array un;                    // u-velocity component in r-direction new
    Array vn;                    // v-velocity component in theta-direction new
    Array wn;                    // w-velocity component in phi-direction new
    Array h2on;                    // water vapour new
    Array h2o_cloudn;                    // water vapour new
    Array h2o_icen;                // cloud water new
    Array h2sn;                    // water vapour new
    Array nh3n;                    // nh3 new
    Array nh3_cloudn;            // nh3_cloud new
    Array nh3_icen;                 // nh3_ice new
    Array ch4n;                    // ch4 new
    Array ch4_cloudn;            // ch4_cloud new
    Array ch4_icen;                 // ch4_ice new
    Array nh4shn;                    // nh4sh new

    Array massflux_h2s;   // mass flux h2s
    Array massflux_nh3;   // mass flux nh3
    Array massflux_nh4sh;   // mass flux nh4sh

    Array difflux_h2s;   // diffusive flux h2s
    Array difflux_nh3;   // diffusive flux nh3
    Array difflux_nh4sh;   // diffusive flux nh4sh

    Array thermalmassflux;   // thermal massflux_h2s

    // Hydrostatic pressure perturbation, the r-integral of the buoyancy. Filled once per
    // Runge-Kutta step by computeHydrostaticPressure(); see the note there.
    Array p_hydro;

    Array p_dyn;                // dynamic pressure
    Array p_dynn;                // dynamic pressure
    Array p_stat;                // static pressure

    Array rhs_t;                // auxilliar field RHS temperature
    Array rhs_u;                // auxilliar field RHS u-velocity component
    Array rhs_v;                // auxilliar field RHS v-velocity component
    Array rhs_w;                // auxilliar field RHS w-velocity component

    Array rhs_h2o;                // auxilliar field RHS water vapour
    Array rhs_h2o_cloud;            // auxilliar field RHS cloud water
    Array rhs_h2o_ice;                // auxilliar field RHS cloud ice
    Array rhs_h2s;                // auxilliar field RHS water vapour
    Array rhs_h2s_cloud;            // auxilliar field RHS cloud water
    Array rhs_h2s_ice;                // auxilliar field RHS cloud ice
    Array rhs_nh3;                // auxilliar field RHS nh3
    Array rhs_nh3_cloud;        // auxilliar field RHS nh3_cloud
    Array rhs_nh3_ice;            // auxilliar field RHS nh3_ice
    Array rhs_ch4;                // auxilliar field RHS ch4
    Array rhs_ch4_cloud;        // auxilliar field RHS ch4_cloud
    Array rhs_ch4_ice;            // auxilliar field RHS ch4_ice
    Array rhs_nh4sh;                // auxilliar field RHS nh4sh
    Array rhs_tke;                // auxilliar field RHS turbulent kinetic energy k*
    Array rhs_dis;                // auxilliar field RHS dissipation (epsilon* or omega*)

    Array fluxlim_nh4sh;  // TVD flux-limiter correction for nh4sh advection

    Array aux;                // auxilliar field u-velocity component
    Array aux_u;                // auxilliar field u-velocity component
    Array aux_v;                // auxilliar field v-velocity component
    Array aux_w;                // auxilliar field w-velocity component

    Array Q_Latent;                // latent heat
    Array Q_Sensible;            // sensible heat
    Array Q_rad;                // radiative heating rate [W/m3] (RadiationJup)
    Array radiation;            // layer-centre net radiative flux [W/m2] (RadiationJup)
    Array epsilon;                // layer emissivity (RadiationJup)
    Array P_rain;               // H2O rain    precipitation flux [kg/m2/s] (PrecipitationJup)
    Array P_snow;               // H2O snow    precipitation flux [kg/m2/s]
    Array P_graupel;            // H2O graupel precipitation flux [kg/m2/s]
    Array P_nh3_rain;           // NH3 rain    precipitation flux [kg/m2/s]
    Array P_nh3_snow;           // NH3 snow    precipitation flux [kg/m2/s]
    Array P_nh3_graupel;        // NH3 graupel precipitation flux [kg/m2/s]
    Array P_ch4_rain;           // CH4 rain    precipitation flux [kg/m2/s]
    Array P_ch4_snow;           // CH4 snow    precipitation flux [kg/m2/s]
    Array P_ch4_graupel;        // CH4 graupel precipitation flux [kg/m2/s]
    Array P_nh4sh;              // NH4SH crystal sedimentation flux [kg/m2/s]
    Array Q_precip;             // latent heating rate from precip phase changes [W/m3] (diagnostic)

    // Precipitation SOURCE TERMS for the moisture equations, as NONDIMENSIONAL tendencies
    // (physical rate [kg/m3/s] already multiplied by L_atm[m]/u_0, the model's time unit), so
    // RHSJup adds them straight into rhs_*. Filled by PrecipitationJup only when
    // ATJUP_PRECIP_COUPLING is on; zero otherwise. See the note at the writeback in
    // PrecipitationJup.h for why the alternative — writing the depleted field in place — does
    // not survive the Runge-Kutta.
    // NH4SH needs no entry here: its Stokes settling is already a term in rhs_nh4sh, and
    // PrecipitationJup::sedimentNH4SH only diagnoses the flux.
    Array S_precip_h2o;         // vapour source (rain evaporating back)   [nondim tendency]
    Array S_precip_h2o_cloud;   // cloud sink (autoconversion, accretion, riming)
    Array S_precip_h2o_ice;     // ice sink (ice autoconversion)
    Array S_precip_nh3;
    Array S_precip_nh3_cloud;
    Array S_precip_nh3_ice;
    Array S_precip_ch4;
    Array S_precip_ch4_cloud;
    Array S_precip_ch4_ice;

    // --- Turbulence closure (TurbulenceJup.h): k-epsilon / k-omega / k-omega SST ---
    // Mirrors the array set of ATOM's TurbulenceAtm.h, and the same DIMENSIONLESS
    // convention: k* = k/u_0^2, omega* = omega_phys*L_atm/u_0, eps* = eps_phys*L_atm/u_0^3,
    // nue* = nue/(u_0*L_atm). NOTE ATJUP's L_atm is in KILOMETRES where ATOM's is in metres,
    // so TurbulenceJup converts with L_atm*1e3 throughout — see its units note.
    Array tke;                  // turbulent kinetic energy k*
    Array tken;                 // k* at time level n
    Array dis;                  // dissipation: epsilon* (k-eps) or omega* (k-omega, SST)
    Array disn;                 // dis at time level n
    Array nue;                  // turbulent (eddy) viscosity nue*
    Array prod;                 // production tensor contraction P_k
    Array tke_source;           // k  source: production - destruction
    Array dis_source;           // dis source: production - destruction + cross-diffusion
    Array_2D vel_star;          // per-column friction velocity u_tau [m/s]

    double re_turb = 1.0;                       // = vel_star_ref*z_0/nue_air, set by TurbulenceJup
    // turb_model ("none" | "k_epsilon" | "k_omega" | "k_omega_SST") is declared by
    // JupiterParams.h.inc from param.py, defaulting to k_omega_SST; ATJUP_TURB_MODEL overrides it.
    // Boundary-layer depth used by the ABL seeding profile and the eddy-viscosity taper.
    // ATOM's value is ~1500 m of terrestrial ABL. There is no Jovian surface boundary layer,
    // so on ATJUP this only has meaning as the depth of the shear layer above the SeaMount;
    // it is metres, and deliberately decoupled from the grid scale (as in ATOM).
    double abl_height = 20000.0;                // [m]
    Array CoriolisForce;        // Coriolis force
    Array CentrifugalForce;             // centrifugal force
    Array BuoyancyForce;        // buoyancy force, Boussinesque approximation
    Array PresGradForce;// pressure gradient force
    Array SeaMount;             // sea mount contour

    // Wall-adjacent eddy viscosity for the momentum equations, in the same units as 1/re, filled
    // once by computeWallViscosity() after the SeaMount contour exists. Zero everywhere except
    // within a few cells of the obstacle. See the long note on the function in InitValues_Jup.cpp.
    Array wall_nue;

    Array w_nh3;                // reaction rate nh3
    Array w_h2s;                // reaction rate h2s
    Array w_nh4sh;                // reaction rate nh4sh

    Array j_nh3;                // ordinary-diffusion mass flux of nh3
    Array j_h2s;                // ordinary-diffusion mass flux of h2s
    Array j_nh4sh;              // ordinary-diffusion mass flux of nh4sh

    Array jT_nh3;               // thermo-diffusion mass flux of nh3
    Array jT_h2s;               // thermo-diffusion mass flux of h2s
    Array jT_nh4sh;             // thermo-diffusion mass flux of nh4sh
};

#endif
