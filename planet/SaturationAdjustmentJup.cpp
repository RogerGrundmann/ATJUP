#include "SaturationAdjustmentJup.h"
#include "cJupiterModel.h"

using namespace std;

void SaturationAdjustmentJup::run(
        const std::string& gas,
        double coeff_A,   double coeff_B,
        double coeff_A_i, double coeff_B_i,
        double t_0,       double t_00,
        double ep,        double lv,  double ls,
        double cp_gas,    double r_gas,
        double C,         double L0,  double R,
        double del_alf,   double del_bet,   double m_mol,
        double C_i,       double L0_i,
        double del_alf_i, double del_bet_i,
        Array& c,         Array& cloud,   Array& ice)
{
    cout << endl << "      SaturationAdjustment of " << gas << " begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // Precompute constant exponent used in the pressure update
    const double exp_pressure = m.g / (m.gam * m.R_ref);
    const double t_range_inv  = 1.0 / (t_0 - t_00);

    // Shared diagnostic state — protected by omp critical where written
    bool   sat_found       = false;
    int    iter_prec_found = 0;
    int    i_sat = 0, j_sat = 0, k_sat = 0;
    double height_sat = 0.0;
    double t_latent = 0.0, p_latent = 0.0;
    double t_sat    = 0.0, p_sat    = 0.0;
    double t_u_sat  = 0.0, p_u_sat  = 0.0;
    double q_v_b_sat = 0.0, q_c_b_sat = 0.0, q_i_b_sat = 0.0;
    double saturation = 0.0;


    // -----------------------------------------------------------------------
    // Main saturation-adjustment loop — fully independent per cell
    // -----------------------------------------------------------------------
    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < m.km; k++){
        for(int j = 0; j < m.jm; j++){
            for(int i = 0; i < m.im; i++){

                const double t_u = m.t.x[i][j][k] * m.t_ref;
                const double p_u = m.p_stat.x[i][j][k];

                // Skip cells that carry no thermodynamic state: the solid SeaMount interior,
                // where bcSolidGround sets t = p_stat = 0, and any cell that has lost a
                // positive temperature or pressure.
                //
                // Both formulas below are singular there. saturation_vapour_pressure() forms
                // -L0/T_K (= -inf at T=0) and del_alf*log(T_K) (= 0*-inf = NaN), and q_Rain_0
                // divides by p_u. The existing `q_Rain_0 <= 0.0` skip cannot catch that: every
                // comparison against NaN is false, so both that test and `c <= q_Rain_0` fail
                // and the cell ran the full mixed-phase iteration on NaN, writing NaN back into
                // c/cloud/ice — one of the two seeds of the domain-wide NaN (the other was the
                // Magnus/Clausius-Clapeyron mix-up in Thermo_Jup.cpp::Latent_Heat).
                //
                // The test is written as !(x > 0.0) rather than (x <= 0.0) so that a NaN that
                // has already arrived from elsewhere is skipped too instead of being propagated.
                if(!(t_u > 0.0) || !(p_u > 0.0) || m.SeaMount.x[i][j][k] == 1.0) continue;

                // Enforce physical bounds
                if(t_u > t_0)              ice.x[i][j][k]   = 0.0;
                if(c.x[i][j][k]     < 0.0) c.x[i][j][k]     = 0.0;
                if(cloud.x[i][j][k] < 0.0) cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k]   < 0.0) ice.x[i][j][k]   = 0.0;

                const double E_Rain_0 = saturation_vapour_pressure(t_u, C, L0, R, del_alf, del_bet);
                const double q_Rain_0 = m.r_mix * ep * E_Rain_0 / p_u;

                // skip: subsaturated, already at saturation, or SVP underflowed to 0 at very cold cells
                if(c.x[i][j][k] <= q_Rain_0 || q_Rain_0 <= 0.0) continue;

                // ---- Mixed-phase iteration (Tao et al. 1988) ----
                double q_v_b   = c.x[i][j][k];
                double q_c_b   = cloud.x[i][j][k];
                double q_i_b   = ice.x[i][j][k];
                double T       = t_u;
                double q_v_hyp = q_v_b;

                bool cell_found = false;
                int  cell_iter  = 0;

                for(int itr = 1; itr <= iter_prec_end; itr++){
                    double CND = (T - t_00) * t_range_inv;
                    double DEP = (t_0 - T)  * t_range_inv;
                    if(T <= t_00){ CND = 0.0; DEP = 1.0; }
                    if(T >= t_0) { CND = 1.0; DEP = 0.0; }

                    const double d_q_v = q_v_hyp - q_v_b;
                    const double d_q_c = -d_q_v * CND;
                    const double d_q_i = -d_q_v * DEP;

                    T     += (lv * d_q_c + ls * d_q_i) / (m.cp_mix * m.r_mix);
                    q_v_b += d_q_v;
                    q_c_b += d_q_c;
                    q_i_b += d_q_i;

                    if(q_v_b < 0.0) q_v_b = 0.0;
                    if(q_c_b < 0.0) q_c_b = 0.0;
                    if(q_i_b < 0.0) q_i_b = 0.0;

                    const double E_Rain = saturation_vapour_pressure(T, C,   L0,   R, del_alf,   del_bet);
                    const double E_Ice  = saturation_vapour_pressure(T, C_i, L0_i, R, del_alf_i, del_bet_i);
                    const double q_Rain = m.r_mix * ep * E_Rain / p_u;
                    const double q_Ice  = m.r_mix * ep * E_Ice  / p_u;

                    if(q_c_b > 0.0 && q_i_b > 0.0)
                        q_v_hyp = (q_c_b * q_Rain + q_i_b * q_Ice) / (q_c_b + q_i_b);
                    else if(q_i_b == 0.0) q_v_hyp = q_Rain;
                    else                  q_v_hyp = q_Ice;

                    if(T >= t_0) q_i_b = 0.0;

                    const double q_diff = std::fabs(q_v_b - q_v_hyp) / (q_v_hyp + 1e-20);
                    if(q_diff <= q_diff_min){
                        cell_found = true;
                        cell_iter  = itr;
                        break;
                    }
                    q_v_hyp = 0.5 * (q_v_hyp + q_v_b);
                }

                // Write converged cell values back
                c.x[i][j][k]     = q_v_b;
                cloud.x[i][j][k] = q_c_b;
                ice.x[i][j][k]   = q_i_b;
                m.t.x[i][j][k]   = T / m.t_ref;

                // Diagnostic capture for the reference reporting cell (equator, mid-atmosphere)
//                if(cell_found && j == m.jm/2 && k == m.km/2 && i == m.im/2){
                if(cell_found){
                    #pragma omp critical
                    {
                        sat_found       = true;
                        iter_prec_found = cell_iter;
                        i_sat = i; j_sat = j; k_sat = k;
                        height_sat = m.get_layer_height(i_sat);
                        t_u_sat    = t_u;
                        p_u_sat    = p_u;
                        t_sat      = T;
                        p_sat      = m.p_ref * std::pow(t_sat / m.t_ref, exp_pressure);
                        t_latent   = t_sat - t_u_sat;
                        p_latent   = p_sat - p_u_sat;
                        q_v_b_sat  = q_v_b;
                        q_c_b_sat  = q_c_b;
                        q_i_b_sat  = q_i_b;
                        saturation = q_v_b - q_Rain_0;
                    }
                }
            }
        }
    }


    // -----------------------------------------------------------------------
    // Clamp negative values left by the iteration
    // -----------------------------------------------------------------------
    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < m.km; k++){
        for(int j = 0; j < m.jm; j++){
            for(int i = 0; i < m.im; i++){
                if(c.x[i][j][k]     <= 0.0) c.x[i][j][k]     = 0.0;
                if(cloud.x[i][j][k] <= 0.0) cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k]   <= 0.0) ice.x[i][j][k]   = 0.0;
            }
        }
    }


    // -----------------------------------------------------------------------
    // Timing and report
    // -----------------------------------------------------------------------
    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for SaturationAdjustment\n",
           elapsed.count() * 1e-9);

    if(!sat_found){
        cout << "      NO saturation found in SaturationAdjustment of " << gas << endl
             << "      iter_prec_end = " << iter_prec_end << endl;
        sat_found = false;
    } else {
        cout << "      saturation of water vapour in SaturationAdjustment of " << gas << " found" << endl
             << "      iter_prec_found = " << iter_prec_found
             << "   iter_prec_end = "      << iter_prec_end   << endl
             << "      i_sat = "     << i_sat
             << "   j_sat = "        << j_sat
             << "   k_sat = "        << k_sat
             << "   height_sat[km] = " << height_sat  << endl
             << "      p_stat[bar] = "  << p_sat
             << "   p_u[bar] = "        << p_u_sat
             << "   p_latent[bar] = "   << p_latent   << endl
             << "      T[°C] = "        << t_sat    - m.t_ref
             << "   t_u[°C] = "         << t_u_sat  - m.t_ref
             << "   t_latent[°C] = "    << t_latent  << endl
             << "      saturation[g/m³] = " << saturation * 1e3 << endl
             << "      " << gas << " humid[g/m³] = "  << q_v_b_sat * 1e3
             << "   cloud[g/m³] = "  << q_c_b_sat * 1e3
             << "   ice[g/m³] = "    << q_i_b_sat * 1e3 << endl;
    }

    cout << "      SaturationAdjustment of " << gas << " ended" << endl;
}
