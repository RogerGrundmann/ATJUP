/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * Thermal-wind initialisation — see the long note in ThermalWindJup.h.
 */

#include "cJupiterModel.h"
#include "ThermalWindJup.h"

using namespace std;

void ThermalWindJup::run(){
    cout << endl << "      ATJUP: ThermalWindJup" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ---- The coefficient that carries drho/dtheta into a nondimensional shear ----
    //
    // The balance is between two terms of rhs_u/rhs_v that carry DIFFERENT nondimensionalisation
    // factors (see the note in RHS_Jup_Turb.cpp): the buoyancy 1e5*L/u_0^2 and the Coriolis L/u_0.
    // Only their ratio survives into the thermal-wind relation, and it is 1e5/u_0 — independent of
    // L, which is why this routine does not need to know the shell thickness at all.
    //
    // B is the model's own buoyancy expression, g*p/(r_mix*R_mix*T*t_ref), read through
    // buoy_pressure() so that it is the same quantity the momentum equation uses, hydrostatic
    // pressure only. B is proportional to rho, so dB/dtheta is drho/dtheta up to that constant.
    const double coeff = 1.0e5 / m.u_0;

    // Latitude band around the equator where geostrophy is abandoned rather than divided by a
    // vanishing f. Expressed as a latitude; |cos(theta)| is the sine of the latitude, so the
    // comparison below is against sin of it.
    static const double lat_taper_deg = [](){
        const char* e = getenv("ATJUP_TW_LAT_TAPER"); return e ? atof(e) : 10.0; }();
    const double cos_taper = std::sin(lat_taper_deg * M_PI / 180.0);

    // Whether to adjust the meridional wind as well, from the ZONAL density gradient. Default
    // OFF, and the measurement is why: it changes v by 15 m/s rms, as much as the zonal
    // adjustment changes w. The initial temperature field is very nearly zonally symmetric, so
    // that is not a broad baroclinic signal — it comes from the obstacle and from whatever
    // longitudinal structure the initialisation leaves behind, and geostrophic balance is the
    // wrong model for the flow around a body. ATJUP_TW_MERIDIONAL=1 enables it.
    static const bool do_v = [](){
        const char* e = getenv("ATJUP_TW_MERIDIONAL"); return e ? atoi(e) != 0 : false; }();

    const double sinthe_min = cJupiterModel::sinthe_min();

    double max_dw = 0.0, max_dv = 0.0, sum_dw2 = 0.0, sum_dv2 = 0.0;
    long long n_cells = 0, n_tapered = 0, n_columns = 0, n_scored = 0;
    double res_before2 = 0.0, res_after2 = 0.0;   // thermal-wind residual, before and after

    std::vector<double> shear_w(m.im), shear_v(m.im), S(m.im), wt(m.im);

    #pragma omp parallel for collapse(2) schedule(dynamic, 4) \
        firstprivate(shear_w, shear_v, S, wt) \
        reduction(+:sum_dw2, sum_dv2, n_cells, n_tapered, n_columns, n_scored, res_before2, res_after2) \
        reduction(max:max_dw, max_dv)
    for(int j = 1; j < m.jm - 1; j++){
        for(int k = 1; k < m.km - 1; k++){

            const double costhe = std::cos(m.the.z[j]);
            const double f      = 2.0 * m.omega * costhe;

            // Blend weight: 1 well away from the equator, 0 at it, quadratic in between so the
            // transition has no kink.
            double alpha = std::fabs(costhe) / cos_taper;
            if(alpha > 1.0) alpha = 1.0;
            alpha = alpha * alpha;
            if(alpha <= 0.0) continue;
            if(alpha < 1.0) n_tapered++;

            double sinthe = std::sin(m.the.z[j]);
            if(sinthe < sinthe_min) sinthe = sinthe_min;   // the metric's own floor, as elsewhere

            // ---- Fluid part of the column ----
            int i0 = 0;
            while(i0 < m.im && m.SeaMount.x[i0][j][k] == 1.0) i0++;
            int i1 = i0;
            while(i1 + 1 < m.im && m.SeaMount.x[i1 + 1][j][k] != 1.0) i1++;
            if(i1 - i0 < 2) continue;

            // ---- The shear the density field demands, level by level ----
            bool ok = true;
            for(int i = i0; i <= i1; i++){
                auto B = [&](int jj, int kk){
                    const double T = m.t.x[i][jj][kk];
                    if(!(T > 0.0)) return 0.0;
                    return m.g * m.buoy_pressure(i, jj, kk) / (m.r_mix * m.R_mix * T * m.t_ref);
                };
                const double dBdthe = (B(j+1, k) - B(j-1, k)) / (2.0 * m.dthe);
                const double dBdphi = (B(j, k+1) - B(j, k-1)) / (2.0 * m.dphi);
                const double rm = m.rad.z[i];

                // The taper multiplies the SHEAR, not just the final blend. 1/f grows without
                // bound towards the equator, so a shear left untapered is astronomically large
                // there however little of it is blended in — and it poisons any diagnostic built
                // from it. With the factor here the product alpha/f behaves as cos(theta) and
                // goes to zero with the Coriolis parameter, which is the physically honest
                // statement: there is no thermal wind at the equator, not a very large one.
                shear_w[i] =  alpha * coeff * dBdthe / (f * rm);
                shear_v[i] = -alpha * coeff * dBdphi / (f * rm * sinthe);
                if(!std::isfinite(shear_w[i]) || !std::isfinite(shear_v[i])) ok = false;

                // dp weights, as in ConvectiveAdjustmentJup: the barotropic mode of a jet is its
                // MASS-weighted mean, and this shell spans 11 bar to 0.02, so an arithmetic mean
                // would not preserve it.
                const double p_lo = (i > i0) ? 0.5 * (m.p_stat.x[i-1][j][k] + m.p_stat.x[i][j][k])
                                             : m.p_stat.x[i0][j][k];
                const double p_hi = (i < i1) ? 0.5 * (m.p_stat.x[i][j][k] + m.p_stat.x[i+1][j][k])
                                             : m.p_stat.x[i1][j][k];
                wt[i] = p_lo - p_hi;
                if(!(wt[i] > 0.0)) ok = false;
            }
            if(!ok) continue;

            // Residual before: how far the model's own shear is from the demanded one. Counted
            // ONLY outside the equatorial taper (alpha == 1). Inside it the relation is not being
            // imposed and is not meaningful, and including those columns would swamp the average
            // with the 1/f the taper exists to avoid.
            const bool score = (alpha >= 1.0);
            if(score){
                for(int i = i0 + 1; i < i1; i++){
                    const double dwdr = (m.w.x[i+1][j][k] - m.w.x[i-1][j][k]) / (2.0 * m.dr);
                    const double r    = dwdr - shear_w[i];
                    res_before2 += r * r;
                    n_scored++;
                }
            }

            // ---- Integrate the shear, then re-centre on the model's own barotropic jet ----
            // Applied to the zonal wind first, then the meridional one with its own integral.
            for(int pass = 0; pass < (do_v ? 2 : 1); pass++){
                const std::vector<double>& sh = (pass == 0) ? shear_w : shear_v;
                Array& fld = (pass == 0) ? m.w : m.v;

                S[i0] = 0.0;
                for(int i = i0 + 1; i <= i1; i++)
                    S[i] = S[i-1] + 0.5 * (sh[i-1] + sh[i]) * m.dr;

                double wsum = 0.0, mean_S = 0.0, mean_f = 0.0;
                for(int i = i0; i <= i1; i++){
                    wsum   += wt[i];
                    mean_S += wt[i] * S[i];
                    mean_f += wt[i] * fld.x[i][j][k];
                }
                mean_S /= wsum;
                mean_f /= wsum;

                for(int i = i0; i <= i1; i++){
                    const double old_v = fld.x[i][j][k];
                    const double bal   = mean_f + (S[i] - mean_S);
                    const double new_v = alpha * bal + (1.0 - alpha) * old_v;
                    const double d     = std::fabs(new_v - old_v) * m.u_0;
                    if(pass == 0){
                        if(d > max_dw) max_dw = d;
                        sum_dw2 += d * d;
                        n_cells++;
                    } else {
                        if(d > max_dv) max_dv = d;
                        sum_dv2 += d * d;
                    }
                    fld.x[i][j][k] = new_v;
                }
            }

            if(score){
                for(int i = i0 + 1; i < i1; i++){
                    const double dwdr = (m.w.x[i+1][j][k] - m.w.x[i-1][j][k]) / (2.0 * m.dr);
                    const double r    = dwdr - shear_w[i];
                    res_after2 += r * r;
                }
            }

            n_columns++;
        }
    }

    const double rms_dw = (n_cells > 0) ? std::sqrt(sum_dw2 / double(n_cells)) : 0.0;
    const double rms_dv = (n_cells > 0) ? std::sqrt(sum_dv2 / double(n_cells)) : 0.0;
    // The residual is a nondimensional shear; multiplied by u_0 it is m/s per unit of r, i.e. per
    // shell thickness, which is the form the measurement in the header is quoted in.
    const double res_b = (n_scored > 0) ? std::sqrt(res_before2 / double(n_scored)) * m.u_0 : 0.0;
    const double res_a = (n_scored > 0) ? std::sqrt(res_after2  / double(n_scored)) * m.u_0 : 0.0;

    printf("      ATJUP: thermal wind — %lld columns (%lld inside the %.0f deg equatorial taper),"
           " dw rms %.3f max %.3f m/s, dv rms %.3f max %.3f m/s\n",
           n_columns, n_tapered, lat_taper_deg, rms_dw, max_dw, rms_dv, max_dv);
    printf("      ATJUP: thermal-wind residual rms %.4f -> %.4f m/s per shell (%.1f %% removed)\n",
           res_b, res_a, (res_b > 0.0) ? 100.0 * (1.0 - res_a / res_b) : 0.0);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for ThermalWindJup\n", elapsed.count() * 1e-9);
    cout << "      ATJUP: ThermalWindJup ended" << endl;
}
