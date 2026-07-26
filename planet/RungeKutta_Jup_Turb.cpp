/*
 * Jupiter Atmosphere Circulation Model (ATJUP)
 * 4th order Runge-Kutta integration of the prognostic fields assembled in RHS_Jup_Turb.cpp
 *
 * Renamed from RungeKutta_Jup.cpp when the turbulent kinetic energy k* and its dissipation
 * dis* (epsilon* for k-epsilon, omega* for k-omega / SST) joined the integrated set, mirroring
 * ATOM_Precipitation/atmosphere/RungeKutta_Atm_Turb.cpp. Previously TurbulenceJup::advance()
 * integrated them on its own with a Patankar splitting and without any transport term; they are
 * now two more scalars of this RK4 system, so they see the same advection, the same turbulent
 * diffusion and the same four-stage time integration as temperature and the species.
*/

#include "cJupiterModel.h"

#include <cstdint>
#include <cstring>

using namespace std;

void cJupiterModel::RungeKuttaJup(){
    cout << endl << "      ATJUP: RungeKuttaJup" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ---- k* ceiling, as in ATOM's turbulent RK4 ----
    // The k production term is linear in k while its sink beta*.k.omega is linear too, so an
    // unbalanced production region can grow k exponentially; ATOM hit 1e98 m2/s2 at a single
    // polar cell before capping. The cap is a runaway guard, deliberately far above anything
    // physical: the closure's own Jovian equilibrium is k ~ 67 m2/s2 (see ParaView_Jup.cpp),
    // and 1000 m2/s2 is above the strongest hurricane-core TKE on record. Override with
    // ATJUP_TKE_MAX [m2/s2] if a storm study needs more headroom.
    static const double tke_max_phys = [](){
        const char* e = getenv("ATJUP_TKE_MAX"); return e ? atof(e) : 1000.0; }();
    const double tke_max_nd = tke_max_phys / (u_0 * u_0);
    constexpr double dis_min_nd = 1.0e-10;    // matches TurbulenceJup::dis_min

    // NaN-safe clamp. std::min/std::max compare with unguarded `<`; under -ffast-math
    // (-ffinite-math-only) the compiler may assume the operands are finite, so a NaN can slip
    // through a plain clamp. Detect non-finite values from the IEEE-754 exponent bits instead
    // and fall back to the lower bound. Same trick ATOM uses in RungeKutta_Atm_Turb.cpp.
    auto safe_clamp = [](double v, double lo, double hi) -> double {
        std::uint64_t bits;
        std::memcpy(&bits, &v, sizeof(bits));
        if((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) return lo;
        return (v < lo) ? lo : ((v > hi) ? hi : v);
    };

    // Precompute sin/cos tables — depend only on j
    // sinthe is clamped to a minimum to prevent 1/sin²θ blow-up near the poles.
    // The sequential k-loop in the RK creates an asymmetric phi Laplacian whose
    // error is amplified by inv_rm2sinthe2; clamping keeps the amplification < 1.
    // PressureSolverJup uses the same threshold. Raised 0.4 -> 0.55 (ATOM parity, metric
    // floor ~57°) to curb the polar 1/sin²θ amplification that seeded a long-run pole blow-up.
    constexpr double sinthe_min = 0.55;
    std::vector<double> sinthe_tbl(jm), costhe_tbl(jm);
    for(int j = 0; j < jm; j++){
        sinthe_tbl[j] = std::max(sinthe_min, std::abs(sin(the.z[j])));
        costhe_tbl[j] = cos(the.z[j]);
    }

    const double inv_2dr   = 1.0 / (2.0 * dr);
    const double inv_2dthe = 1.0 / (2.0 * dthe);
    const double inv_2dphi = 1.0 / (2.0 * dphi);
    const double inv_dr2   = 1.0 / (dr   * dr);
    const double inv_dthe2 = 1.0 / (dthe * dthe);
    const double inv_dphi2 = 1.0 / (dphi * dphi);

    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < im-1; i++){
        for(int j = 3; j < jm-3; j++){

            // Build geometry struct once per (i,j)
            CellGeometry geo;
            geo.rm      = rad.z[i];
            geo.rm2     = geo.rm * geo.rm;
            geo.exp_rm   = coord_stretching ? 1.0 / (geo.rm + 1.0) : 1.0;
            geo.exp_2_rm = geo.exp_rm * geo.exp_rm;
            geo.sinthe  = sinthe_tbl[j];
            geo.sinthe2 = geo.sinthe * geo.sinthe;
            geo.costhe  = costhe_tbl[j];
            if(j > 90) geo.costhe = -geo.costhe;
            geo.cotanthe            = geo.costhe / geo.sinthe;
            geo.inv_rm              = 1.0 / geo.rm;
            geo.inv_rm2             = 1.0 / geo.rm2;
            geo.inv_rmsinthe        = 1.0 / (geo.rm * geo.sinthe);
            geo.inv_rm2sinthe       = geo.inv_rm2 / geo.sinthe;
            geo.inv_rm2sinthe2      = geo.inv_rm2 / geo.sinthe2;
            geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;
            geo.inv_2dr   = inv_2dr;
            geo.inv_2dthe = inv_2dthe;
            geo.inv_2dphi = inv_2dphi;
            geo.inv_dr2   = inv_dr2;
            geo.inv_dthe2 = inv_dthe2;
            geo.inv_dphi2 = inv_dphi2;

            for(int k = 1; k < km-1; k++){

                if(SeaMount.x[i][j][k] == 1.0) continue;

                double tn_ijk     = tn.x[i][j][k];
                double un_ijk     = un.x[i][j][k];
                double vn_ijk     = vn.x[i][j][k];
                double wn_ijk     = wn.x[i][j][k];
                double h2on_ijk   = h2on.x[i][j][k];
                double h2ocn_ijk  = h2o_cloudn.x[i][j][k];
                double h2oin_ijk  = h2o_icen.x[i][j][k];
                double h2sn_ijk   = h2sn.x[i][j][k];
                double nh3n_ijk   = nh3n.x[i][j][k];
                double nh3cn_ijk  = nh3_cloudn.x[i][j][k];
                double nh3in_ijk  = nh3_icen.x[i][j][k];
                double ch4n_ijk   = ch4n.x[i][j][k];
                double ch4cn_ijk  = ch4_cloudn.x[i][j][k];
                double ch4in_ijk  = ch4_icen.x[i][j][k];
                double nh4shn_ijk = nh4shn.x[i][j][k];
                double tken_ijk   = tken.x[i][j][k];
                double disn_ijk   = disn.x[i][j][k];

                // ----- RK stage 1 -----
                cJupiterModel::RHSJup(i, j, k, geo);
                double kt1     = rhs_t.x[i][j][k];
                double ku1     = rhs_u.x[i][j][k];
                double kv1     = rhs_v.x[i][j][k];
                double kw1     = rhs_w.x[i][j][k];
                double kc1     = rhs_h2o.x[i][j][k];
                double kcloud1 = rhs_h2o_cloud.x[i][j][k];
                double kice1   = rhs_h2o_ice.x[i][j][k];
                double kh2s1   = rhs_h2s.x[i][j][k];
                double knh31   = rhs_nh3.x[i][j][k];
                double knh3c1  = rhs_nh3_cloud.x[i][j][k];
                double knh3i1  = rhs_nh3_ice.x[i][j][k];
                double kch41   = rhs_ch4.x[i][j][k];
                double kch4c1  = rhs_ch4_cloud.x[i][j][k];
                double kch4i1  = rhs_ch4_ice.x[i][j][k];
                double knh4sh1 = rhs_nh4sh.x[i][j][k];
                double ktke1   = rhs_tke.x[i][j][k];
                double kdis1   = rhs_dis.x[i][j][k];

                t.x[i][j][k]         = tn_ijk    + kt1     * 0.5 * dt;
                u.x[i][j][k]         = un_ijk    + ku1     * 0.5 * dt;
                v.x[i][j][k]         = vn_ijk    + kv1     * 0.5 * dt;
                w.x[i][j][k]         = wn_ijk    + kw1     * 0.5 * dt;
                h2o.x[i][j][k]       = h2on_ijk  + kc1     * 0.5 * dt;
                h2o_cloud.x[i][j][k] = h2ocn_ijk + kcloud1 * 0.5 * dt;
                h2o_ice.x[i][j][k]   = h2oin_ijk + kice1   * 0.5 * dt;
                h2s.x[i][j][k]       = h2sn_ijk  + kh2s1   * 0.5 * dt;
                nh3.x[i][j][k]       = nh3n_ijk  + knh31   * 0.5 * dt;
                nh3_cloud.x[i][j][k] = nh3cn_ijk + knh3c1  * 0.5 * dt;
                nh3_ice.x[i][j][k]   = nh3in_ijk + knh3i1  * 0.5 * dt;
                ch4.x[i][j][k]       = ch4n_ijk  + kch41   * 0.5 * dt;
                ch4_cloud.x[i][j][k] = ch4cn_ijk + kch4c1  * 0.5 * dt;
                ch4_ice.x[i][j][k]   = ch4in_ijk + kch4i1  * 0.5 * dt;
                nh4sh.x[i][j][k]     = nh4shn_ijk + knh4sh1 * 0.5 * dt;
                // k* is positive-definite and dis* strictly positive; clamp at every stage so a
                // stiff intermediate can never feed a negative k or a zero omega back into the
                // next RHS evaluation (nue* = k/omega would then be non-finite).
                tke.x[i][j][k]       = safe_clamp(tken_ijk + ktke1 * 0.5 * dt, 0.0, tke_max_nd);
                dis.x[i][j][k]       = std::max(dis_min_nd, disn_ijk + kdis1 * 0.5 * dt);

                // ----- RK stage 2 -----
                cJupiterModel::RHSJup(i, j, k, geo);
                double kt2     = rhs_t.x[i][j][k];
                double ku2     = rhs_u.x[i][j][k];
                double kv2     = rhs_v.x[i][j][k];
                double kw2     = rhs_w.x[i][j][k];
                double kc2     = rhs_h2o.x[i][j][k];
                double kcloud2 = rhs_h2o_cloud.x[i][j][k];
                double kice2   = rhs_h2o_ice.x[i][j][k];
                double kh2s2   = rhs_h2s.x[i][j][k];
                double knh32   = rhs_nh3.x[i][j][k];
                double knh3c2  = rhs_nh3_cloud.x[i][j][k];
                double knh3i2  = rhs_nh3_ice.x[i][j][k];
                double kch42   = rhs_ch4.x[i][j][k];
                double kch4c2  = rhs_ch4_cloud.x[i][j][k];
                double kch4i2  = rhs_ch4_ice.x[i][j][k];
                double knh4sh2 = rhs_nh4sh.x[i][j][k];
                double ktke2   = rhs_tke.x[i][j][k];
                double kdis2   = rhs_dis.x[i][j][k];

                t.x[i][j][k]         = tn_ijk    + kt2     * 0.5 * dt;
                u.x[i][j][k]         = un_ijk    + ku2     * 0.5 * dt;
                v.x[i][j][k]         = vn_ijk    + kv2     * 0.5 * dt;
                w.x[i][j][k]         = wn_ijk    + kw2     * 0.5 * dt;
                h2o.x[i][j][k]       = h2on_ijk  + kc2     * 0.5 * dt;
                h2o_cloud.x[i][j][k] = h2ocn_ijk + kcloud2 * 0.5 * dt;
                h2o_ice.x[i][j][k]   = h2oin_ijk + kice2   * 0.5 * dt;
                h2s.x[i][j][k]       = h2sn_ijk  + kh2s2   * 0.5 * dt;
                nh3.x[i][j][k]       = nh3n_ijk  + knh32   * 0.5 * dt;
                nh3_cloud.x[i][j][k] = nh3cn_ijk + knh3c2  * 0.5 * dt;
                nh3_ice.x[i][j][k]   = nh3in_ijk + knh3i2  * 0.5 * dt;
                ch4.x[i][j][k]       = ch4n_ijk  + kch42   * 0.5 * dt;
                ch4_cloud.x[i][j][k] = ch4cn_ijk + kch4c2  * 0.5 * dt;
                ch4_ice.x[i][j][k]   = ch4in_ijk + kch4i2  * 0.5 * dt;
                nh4sh.x[i][j][k]     = nh4shn_ijk + knh4sh2 * 0.5 * dt;
                tke.x[i][j][k]       = safe_clamp(tken_ijk + ktke2 * 0.5 * dt, 0.0, tke_max_nd);
                dis.x[i][j][k]       = std::max(dis_min_nd, disn_ijk + kdis2 * 0.5 * dt);

                // ----- RK stage 3 -----
                cJupiterModel::RHSJup(i, j, k, geo);
                double kt3     = rhs_t.x[i][j][k];
                double ku3     = rhs_u.x[i][j][k];
                double kv3     = rhs_v.x[i][j][k];
                double kw3     = rhs_w.x[i][j][k];
                double kc3     = rhs_h2o.x[i][j][k];
                double kcloud3 = rhs_h2o_cloud.x[i][j][k];
                double kice3   = rhs_h2o_ice.x[i][j][k];
                double kh2s3   = rhs_h2s.x[i][j][k];
                double knh33   = rhs_nh3.x[i][j][k];
                double knh3c3  = rhs_nh3_cloud.x[i][j][k];
                double knh3i3  = rhs_nh3_ice.x[i][j][k];
                double kch43   = rhs_ch4.x[i][j][k];
                double kch4c3  = rhs_ch4_cloud.x[i][j][k];
                double kch4i3  = rhs_ch4_ice.x[i][j][k];
                double knh4sh3 = rhs_nh4sh.x[i][j][k];
                double ktke3   = rhs_tke.x[i][j][k];
                double kdis3   = rhs_dis.x[i][j][k];

                t.x[i][j][k]         = tn_ijk    + kt3     * dt;
                u.x[i][j][k]         = un_ijk    + ku3     * dt;
                v.x[i][j][k]         = vn_ijk    + kv3     * dt;
                w.x[i][j][k]         = wn_ijk    + kw3     * dt;
                h2o.x[i][j][k]       = h2on_ijk  + kc3     * dt;
                h2o_cloud.x[i][j][k] = h2ocn_ijk + kcloud3 * dt;
                h2o_ice.x[i][j][k]   = h2oin_ijk + kice3   * dt;
                h2s.x[i][j][k]       = h2sn_ijk  + kh2s3   * dt;
                nh3.x[i][j][k]       = nh3n_ijk  + knh33   * dt;
                nh3_cloud.x[i][j][k] = nh3cn_ijk + knh3c3  * dt;
                nh3_ice.x[i][j][k]   = nh3in_ijk + knh3i3  * dt;
                ch4.x[i][j][k]       = ch4n_ijk  + kch43   * dt;
                ch4_cloud.x[i][j][k] = ch4cn_ijk + kch4c3  * dt;
                ch4_ice.x[i][j][k]   = ch4in_ijk + kch4i3  * dt;
                nh4sh.x[i][j][k]     = nh4shn_ijk + knh4sh3 * dt;
                tke.x[i][j][k]       = safe_clamp(tken_ijk + ktke3 * dt, 0.0, tke_max_nd);
                dis.x[i][j][k]       = std::max(dis_min_nd, disn_ijk + kdis3 * dt);

                // ----- RK stage 4 -----
                cJupiterModel::RHSJup(i, j, k, geo);
                double kt4     = rhs_t.x[i][j][k];
                double ku4     = rhs_u.x[i][j][k];
                double kv4     = rhs_v.x[i][j][k];
                double kw4     = rhs_w.x[i][j][k];
                double kc4     = rhs_h2o.x[i][j][k];
                double kcloud4 = rhs_h2o_cloud.x[i][j][k];
                double kice4   = rhs_h2o_ice.x[i][j][k];
                double kh2s4   = rhs_h2s.x[i][j][k];
                double knh34   = rhs_nh3.x[i][j][k];
                double knh3c4  = rhs_nh3_cloud.x[i][j][k];
                double knh3i4  = rhs_nh3_ice.x[i][j][k];
                double kch44   = rhs_ch4.x[i][j][k];
                double kch4c4  = rhs_ch4_cloud.x[i][j][k];
                double kch4i4  = rhs_ch4_ice.x[i][j][k];
                double knh4sh4 = rhs_nh4sh.x[i][j][k];
                double ktke4   = rhs_tke.x[i][j][k];
                double kdis4   = rhs_dis.x[i][j][k];

                // ----- Final RK4 update -----
                const double one_sixth = 1.0 / 6.0;
                t.x[i][j][k]         = tn_ijk    + dt * (kt1     + 2.0*kt2     + 2.0*kt3     + kt4    ) * one_sixth;
                u.x[i][j][k]         = un_ijk    + dt * (ku1     + 2.0*ku2     + 2.0*ku3     + ku4    ) * one_sixth;
                v.x[i][j][k]         = vn_ijk    + dt * (kv1     + 2.0*kv2     + 2.0*kv3     + kv4    ) * one_sixth;
                w.x[i][j][k]         = wn_ijk    + dt * (kw1     + 2.0*kw2     + 2.0*kw3     + kw4    ) * one_sixth;
                h2o.x[i][j][k]       = std::max(0.0, h2on_ijk  + dt * (kc1     + 2.0*kc2     + 2.0*kc3     + kc4    ) * one_sixth);
                h2o_cloud.x[i][j][k] = std::max(0.0, h2ocn_ijk + dt * (kcloud1 + 2.0*kcloud2 + 2.0*kcloud3 + kcloud4) * one_sixth);
                h2o_ice.x[i][j][k]   = std::max(0.0, h2oin_ijk + dt * (kice1   + 2.0*kice2   + 2.0*kice3   + kice4  ) * one_sixth);
                h2s.x[i][j][k]       = std::max(0.0, h2sn_ijk  + dt * (kh2s1   + 2.0*kh2s2   + 2.0*kh2s3   + kh2s4  ) * one_sixth);
                nh3.x[i][j][k]       = std::max(0.0, nh3n_ijk  + dt * (knh31   + 2.0*knh32   + 2.0*knh33   + knh34  ) * one_sixth);
                nh3_cloud.x[i][j][k] = std::max(0.0, nh3cn_ijk + dt * (knh3c1  + 2.0*knh3c2  + 2.0*knh3c3  + knh3c4 ) * one_sixth);
                nh3_ice.x[i][j][k]   = std::max(0.0, nh3in_ijk + dt * (knh3i1  + 2.0*knh3i2  + 2.0*knh3i3  + knh3i4 ) * one_sixth);
                ch4.x[i][j][k]       = std::max(0.0, ch4n_ijk  + dt * (kch41   + 2.0*kch42   + 2.0*kch43   + kch44  ) * one_sixth);
                ch4_cloud.x[i][j][k] = std::max(0.0, ch4cn_ijk + dt * (kch4c1  + 2.0*kch4c2  + 2.0*kch4c3  + kch4c4 ) * one_sixth);
                ch4_ice.x[i][j][k]   = std::max(0.0, ch4in_ijk + dt * (kch4i1  + 2.0*kch4i2  + 2.0*kch4i3  + kch4i4 ) * one_sixth);
                nh4sh.x[i][j][k]     = std::max(0.0, nh4shn_ijk + dt * (knh4sh1 + 2.0*knh4sh2 + 2.0*knh4sh3 + knh4sh4) * one_sixth);
                tke.x[i][j][k]       = safe_clamp(tken_ijk + dt * (ktke1 + 2.0*ktke2 + 2.0*ktke3 + ktke4) * one_sixth,
                                                  0.0, tke_max_nd);
                dis.x[i][j][k]       = std::max(dis_min_nd,
                                                disn_ijk + dt * (kdis1 + 2.0*kdis2 + 2.0*kdis3 + kdis4) * one_sixth);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for solveRungeKutta\n", elapsed.count() * 1e-9);

    cout << "      ATJUP: RungeKuttaJup ended" << endl;
    return;
}
/*
*
*/
