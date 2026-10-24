/*
 * ATJUP_CHROM=1 (2026-10-24, default 0) -- the CHROMOPHORE TRACER: what colours Jupiter.
 *
 * Every cloud this model condenses (NH3 ice, NH4SH, H2O) is white. The reds, browns and oranges of the
 * photographs come from trace amounts of coloured compounds, the chromophores, made where solar UV breaks
 * up NH3 (with PH3, H2S, C2H2) high in the atmosphere and the fragments recombine. Their identity is still
 * argued (sulphur / phosphorus compounds; NH3 + acetylene products, Carlson et al. 2016); what the models
 * agree on is WHERE they form -- in the UV-lit top, from whatever NH3 the dynamics brings up there -- and
 * that they bleach / are destroyed again. The GRS is red because it lifts NH3 high and keeps it there.
 *
 * So this is a passive species chrom [kg/m3], carried and diffused like the other species, with
 *
 *     d chrom/dt |_chem = J * c_NH3 * F_sun(theta) * A_uv  -  chrom / tau
 *
 *     c_NH3    NH3 vapour + cloud + ice (ATJUP_CHROM_SRC=1, default) or vapour alone (0)
 *     F_sun    sin(colatitude): the diurnal-mean insolation shape, no season (1 at the equator)
 *     A_uv     exp(-N_above / N_uv): the UV left after the NH3 VAPOUR column N_above [kg/m2] above the
 *              cell (from the cell's upper face to the lid) has absorbed it -- NH3 shields itself.
 *              N_uv = ATJUP_CHROM_NUV [kg/m2], default 3e-4: tau = 1 at ~1e22 molecules/m2, i.e. the
 *              NH3 cross-section ~1e-18 cm2 at 190-220 nm. On the TV26 state that puts the UV edge at
 *              i ~27-28, the top of the NH3 ice deck, and deeper over the vortex (less vapour above it).
 *              A fixed decay length (the first version, H_uv 10 km) let the NH3 that grows with depth
 *              win: the tracer peaked at i 26, under the cloud top.
 *     J        ATJUP_CHROM_J, effective production rate [1/s], default 1e-6
 *     tau      ATJUP_CHROM_TAU [days of 86400 s], default 30
 *
 * A_uv is built once per iteration (chromUV(), before RungeKuttaJup) into chrom_uv; within the four RK
 * stages NH3 changes far too little to matter.
 *
 * UNCALIBRATED ON PURPOSE: J lumps photolysis rate, chemical yield and optical strength, none of which
 * is known well enough to fix an absolute amount; read chrom as a relative colour intensity. It takes
 * nothing from NH3 (the yield is tiny) and feeds nothing back (no radiation, no density). Not modelled:
 * haze / cloud scattering in A_uv, settling with the cloud particles, scavenging by precipitation, a
 * temperature-dependent destruction at depth -- the lifetime and the diffusion are its only sinks.
 *
 * Off: chrom is never touched and stays 0, no output carries it and the restart file is unchanged, so a
 * default run is byte-identical. On: the restart file carries chrom and chromn after the other arrays; a
 * restart file without them starts chrom from 0 (load_state says so).
 */

#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace JupChrom {
    struct Par { double J, N_uv, tau_s; int src; };
    inline const Par& par(){
        static const Par p = [](){
            auto d = [](const char* n, double v){ const char* e = getenv(n); return e ? atof(e) : v; };
            return Par{ d("ATJUP_CHROM_J", 1.0e-6), d("ATJUP_CHROM_NUV", 3.0e-4), d("ATJUP_CHROM_TAU", 30.0) * 86400.0,
                        (int)d("ATJUP_CHROM_SRC", 1.0) };
        }();
        return p;
    }
}

// chrom_uv = exp(-N_above / N_uv), N_above the NH3 vapour column from the upper face of the cell to the lid
// (the lid layer im-1 has no thickness, as in every column sum here).
inline void cJupiterModel::chromUV(){
    if(!chrom_on()) return;
    const double n_uv = JupChrom::par().N_uv;
    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            double n_above = 0.0;
            chrom_uv.x[im-1][j][k] = 1.0;
            for(int i = im - 2; i >= 0; i--){
                chrom_uv.x[i][j][k] = n_uv > 0.0 ? std::exp(-n_above / n_uv) : 1.0;
                n_above += std::max(0.0, nh3.x[i][j][k]) * layer_thickness_m(i);
            }
        }
    }
}

inline double cJupiterModel::chromSource(int i, int j, int k){
    const JupChrom::Par& p = JupChrom::par();
    const double c = nh3.x[i][j][k] + (p.src ? nh3_cloud.x[i][j][k] + nh3_ice.x[i][j][k] : 0.0);
    const double f_sun = std::max(0.0, std::sin(the.z[j]));
    const double rate = p.J * std::max(0.0, c) * f_sun * chrom_uv.x[i][j][k]
                      - (p.tau_s > 0.0 ? chrom.x[i][j][k] / p.tau_s : 0.0);     // kg/m3/s
    return rate * (L_atm * 1.0e3 / u_0);                                          // per nondimensional time
}

inline void cJupiterModel::chromInfo(){
    if(!chrom_on()) return;
    const JupChrom::Par& p = JupChrom::par();
    printf("      ATJUP: ATJUP_CHROM on -- chromophore tracer: J %.3e 1/s, N_uv %.3e kg/m2 (NH3 vapour self-shielding), tau %.1f d, source %s\n",
           p.J, p.N_uv, p.tau_s / 86400.0, p.src ? "NH3 vapour + cloud + ice" : "NH3 vapour");
}
