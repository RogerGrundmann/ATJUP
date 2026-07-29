#pragma once

#include "cJupiterModel.h"
#include "Utils.h"

#include <iostream>
#include <vector>
#include <chrono>
#include <cstdlib>   // getenv/atof for ATJUP_U_INIT_SCALE

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace JupiterUtils;

// ============================================================================
// Header-only friend class: initialises u, v, w for Jupiter's 8-cell
// zonal circulation (4 pairs of Hadley/Ferrel cells per hemisphere).
// Accessed from cJupiterModel::JupiterCellStructure().
// ============================================================================
class VelocityInitializerJup {
public:
    explicit VelocityInitializerJup(cJupiterModel& model)
        : m(model)
    {}

    // ========================================================================
    // Main entry point — mirrors cJupiterModel::JupiterCellStructure()
    // ========================================================================
    void compute()
    {
        std::cout << "\n\n\n      ATJUP: VelocityInitializerJup::compute" << std::endl;
        printf("      ATJUP: ATJUP_U_INIT_SCALE = %.4g  ->  radial cell branch peaks at %.3g m/s\n",
               u_init_scale(), 0.5 * 80.0 * u_init_scale());

        auto begin = std::chrono::high_resolution_clock::now();

        // --------------------------------------------------------------------
        // u-velocity: seed tropopause profiles at cell boundaries
        // --------------------------------------------------------------------
        // northern hemisphere cell boundaries (equator → north pole)
        init_u(m.u,   0);  init_u(m.u,   4);  init_u(m.u,   8);  init_u(m.u,  12);   // polar cell (centre j8)
        init_u(m.u,  15);  // 8. Ferrel cell north, min
        init_u(m.u,  16);  // 8. Ferrel cell north, min
        init_u(m.u,  17);  // 8. Ferrel / 8. Hadley center
        init_u(m.u,  19);  // 7. Hadley cell north, max
        // init_u 21 removed (gap closed)
        // init_u 23 removed
        // init_u(m.u, 24) removed — j24 updraft deleted, interpolated across as a sink
        // init_u 25 removed
        // init_u 27 removed (gap closed)
        init_u(m.u,  29);  // 6. Ferrel cell north, min
        init_u(m.u,  30);  // 6. Ferrel / 5. Hadley center
        init_u(m.u,  34);  // 5. Ferrel cell north, min
        init_u(m.u,  35);  // 5. Ferrel center
        init_u(m.u,  37);  // 5. Hadley cell north, max
        init_u(m.u,  40);  // 5. Ferrel / 4. Hadley center
        init_u(m.u,  43);  // 4. Hadley cell north, max
        init_u(m.u,  46);  // 4. Ferrel center
        init_u(m.u,  48);  // 4. Ferrel cell north, min
        init_u(m.u,  49);  // 4. Ferrel / 3. Hadley center
        init_u(m.u,  51);  // 3. Hadley cell north, max
        init_u(m.u,  52);  // 3. Ferrel center
        init_u(m.u,  53);  // 3. Ferrel cell north, min
        init_u(m.u,  54);  // 3. Ferrel / 2. Hadley center
        init_u(m.u,  57);  // 2. Ferrel cell north, min
        init_u(m.u,  60);  // 2. Ferrel center
        init_u(m.u,  64);  // 2. Hadley cell north, max
        init_u(m.u,  66);  // 2. Hadley / 1. Ferrel center
        init_u(m.u,  69);  // 1. Ferrel cell north, min
        init_u(m.u,  76);  // 1. Ferrel center
        init_u(m.u,  83);  // 1. Hadley cell north, max
        // equator
        init_u(m.u,  90);
        // southern hemisphere cell boundaries (equator → south pole)
        init_u(m.u,  97);
        init_u(m.u, 102);
        init_u(m.u, 107);
        init_u(m.u, 111);  // 1. Hadley / 2. Ferrel center
        init_u(m.u, 115);
        init_u(m.u, 118);
        init_u(m.u, 121);
        init_u(m.u, 123);  // 2. Ferrel / 3. Hadley center
        init_u(m.u, 124);
        init_u(m.u, 125);
        init_u(m.u, 126);
        init_u(m.u, 127);  // 3. Ferrel / 4. Hadley center
        init_u(m.u, 128);
        init_u(m.u, 129);
        init_u(m.u, 130);
        init_u(m.u, 136);  // 4. Ferrel / 5. Hadley center
        init_u(m.u, 138);
        init_u(m.u, 140);
        init_u(m.u, 142);
        init_u(m.u, 143);  // 5. Ferrel / 6. Hadley center
        init_u(m.u, 145);
        init_u(m.u, 146);
        init_u(m.u, 149);
        init_u(m.u, 152);  // 6. Ferrel / 7. Hadley center
        init_u(m.u, 153);
        init_u(m.u, 154);
        init_u(m.u, 155);
        init_u(m.u, 157);  // 7. Ferrel / 8. Hadley center
        init_u(m.u, 158);
        init_u(m.u, 160);
        init_u(m.u, 166);  init_u(m.u, 168);  init_u(m.u, 170);  init_u(m.u, 172);   // polar cells (2, not mirrored)
        init_u(m.u, 174);  init_u(m.u, 176);  init_u(m.u, 178);  init_u(m.u, 180);

        // --------------------------------------------------------------------
        // v/w velocity: linear profile between surface and tropopause values
        // --------------------------------------------------------------------
        // equator
        init_v_or_w(m.v,  90,    0.0,   0.0);
        init_v_or_w(m.w,  90,   73.0,   0.0);

        // 1. Hadley cell
        init_v_or_w(m.v,  83, -40.0,  40.0);   init_v_or_w(m.w,  83, 137.0,   0.0);
        init_v_or_w(m.v,  97,  40.0, -40.0);   init_v_or_w(m.w,  97, 112.0,   0.0);
        // 1. Ferrel cell
        init_v_or_w(m.v,  76,   0.0,   0.0);   init_v_or_w(m.w,  76,  13.0,   0.0);
        init_v_or_w(m.v, 102,   0.0,   0.0);   init_v_or_w(m.w, 102,  40.0,   0.0);
        init_v_or_w(m.v,  69,  40.0, -40.0);   init_v_or_w(m.w,  69, -88.0,   0.0);
        init_v_or_w(m.v, 107, -40.0,  40.0);   init_v_or_w(m.w, 107, -25.0,   0.0);

        // 2. Hadley cell
        init_v_or_w(m.v,  66,   0.0,   0.0);   init_v_or_w(m.w,  66,  -6.0,   0.0);
        init_v_or_w(m.v, 111,   0.0,   0.0);   init_v_or_w(m.w, 111,  40.0,   0.0);
        init_v_or_w(m.v,  64, -40.0,  40.0);   init_v_or_w(m.w,  64,  40.0,   0.0);
        init_v_or_w(m.v, 115,  40.0, -40.0);   init_v_or_w(m.w, 115, 162.0,   0.0);
        // 2. Ferrel cell
        init_v_or_w(m.v,  60,   0.0,   0.0);   init_v_or_w(m.w,  60,  13.0,   0.0);
        init_v_or_w(m.v, 118,   0.0,   0.0);   init_v_or_w(m.w, 118,  13.0,   0.0);
        init_v_or_w(m.v,  57,  40.0, -40.0);   init_v_or_w(m.w,  57, -13.0,   0.0);
        init_v_or_w(m.v, 121, -40.0,  40.0);   init_v_or_w(m.w, 121,  -1.0,   0.0);

        // 3. Hadley cell
        init_v_or_w(m.v,  54,   0.0,   0.0);   init_v_or_w(m.w,  54,   1.0,   0.0);
        init_v_or_w(m.v, 123,   0.0,   0.0);   init_v_or_w(m.w, 123,   1.0,   0.0);
        init_v_or_w(m.v,  53, -40.0,  40.0);   init_v_or_w(m.w,  53,  31.0,   0.0);
        init_v_or_w(m.v, 124,  40.0, -40.0);   init_v_or_w(m.w, 124,  -1.0,   0.0);
        // 3. Ferrel cell
        init_v_or_w(m.v,  52,   0.0,   0.0);   init_v_or_w(m.w,  52,  13.0,   0.0);
        init_v_or_w(m.v, 125,   0.0,   0.0);   init_v_or_w(m.w, 125,  13.0,   0.0);
        init_v_or_w(m.v,  51,  40.0, -40.0);   init_v_or_w(m.w,  51,  -6.0,   0.0);
        init_v_or_w(m.v, 126, -40.0,  40.0);   init_v_or_w(m.w, 126, -12.0,   0.0);

        // 4. Hadley cell
        init_v_or_w(m.v,  49,   0.0,   0.0);   init_v_or_w(m.w,  49,  13.0,   0.0);
        init_v_or_w(m.v, 127,   0.0,   0.0);   init_v_or_w(m.w, 127,   0.0,   0.0);
        init_v_or_w(m.v,  48, -40.0,  40.0);   init_v_or_w(m.w,  48,  33.0,   0.0);
        init_v_or_w(m.v, 128,  40.0, -40.0);   init_v_or_w(m.w, 128,  18.0,   0.0);
        // 4. Ferrel cell
        init_v_or_w(m.v,  46,   0.0,   0.0);   init_v_or_w(m.w,  46,  40.0,   0.0);
        init_v_or_w(m.v, 129,   0.0,   0.0);   init_v_or_w(m.w, 129,  25.0,   0.0);
        init_v_or_w(m.v,  43,  40.0, -40.0);   init_v_or_w(m.w,  43,  -6.0,   0.0);
        init_v_or_w(m.v, 130, -40.0,  40.0);   init_v_or_w(m.w, 130, -12.0,   0.0);

        // 5. Hadley cell
        init_v_or_w(m.v,  40,   0.0,   0.0);   init_v_or_w(m.w,  40,   1.0,   0.0);
        init_v_or_w(m.v, 136,   0.0,   0.0);   init_v_or_w(m.w, 136,   0.0,   0.0);
        init_v_or_w(m.v,  37, -40.0,  40.0);   init_v_or_w(m.w,  37,  37.0,   0.0);
        init_v_or_w(m.v, 138,  40.0, -40.0);   init_v_or_w(m.w, 138,  -6.0,   0.0);
        // 5. Ferrel cell
        init_v_or_w(m.v,  35,   0.0,   0.0);   init_v_or_w(m.w,  35,  13.0,   0.0);
        init_v_or_w(m.v, 140,   0.0,   0.0);   init_v_or_w(m.w, 140,  -6.0,   0.0);
        init_v_or_w(m.v,  34,  40.0, -40.0);   init_v_or_w(m.w,  34,  -6.0,   0.0);
        init_v_or_w(m.v, 142, -40.0,  40.0);   init_v_or_w(m.w, 142, -13.0,   0.0);

        // 6. Hadley cell
        init_v_or_w(m.v,  30,   0.0,   0.0);   init_v_or_w(m.w,  30,   1.0,   0.0);
        init_v_or_w(m.v, 143,   0.0,   0.0);   init_v_or_w(m.w, 143,   1.0,   0.0);
        init_v_or_w(m.v,  29, -40.0,  40.0);   init_v_or_w(m.w,  29,  12.0,   0.0);
        init_v_or_w(m.v, 145,  40.0, -40.0);   init_v_or_w(m.w, 145,  13.0,   0.0);
        // 6. Ferrel cell
        // j27 boundary removed (gap closed: cells j17         init_v_or_w(m.v,  27,   0.0,   0.0);   init_v_or_w(m.w,  27,   1.0,   0.0); j30 merge)
        init_v_or_w(m.v, 146,   0.0,   0.0);   init_v_or_w(m.w, 146,  -6.0,   0.0);
        // j25 flank removed (belonged to deleted j24 cell)
        init_v_or_w(m.v, 149, -40.0,  40.0);   init_v_or_w(m.w, 149,  -1.0,   0.0);

        // 7. Hadley cell
        // j=24 updraft centre deleted (redundant vertical line): neighbours j17 & j30 expand across it
        init_v_or_w(m.v, 152,   0.0,   0.0);   init_v_or_w(m.w, 152,  -1.0,   0.0);
        // j23 flank removed (belonged to deleted j24 cell)
        init_v_or_w(m.v, 153,  40.0, -40.0);   init_v_or_w(m.w, 153,  -1.0,   0.0);
        // 7. Ferrel cell
        // j21 boundary removed (gap closed: cells j17         init_v_or_w(m.v,  21,   0.0,   0.0);   init_v_or_w(m.w,  21,  33.0,   0.0); j30 merge)
        init_v_or_w(m.v, 154,   0.0,   0.0);   init_v_or_w(m.w, 154,  13.0,   0.0);
        init_v_or_w(m.v,  19,  40.0, -40.0);   init_v_or_w(m.w,  19,  -6.0,   0.0);
        init_v_or_w(m.v, 155, -40.0,  40.0);   init_v_or_w(m.w, 155,  -6.0,   0.0);

        // 8. Ferrel cell
        init_v_or_w(m.v,  17,   0.0,   0.0);   init_v_or_w(m.w,  17,   0.0,   0.0);
        init_v_or_w(m.v, 157,   0.0,   0.0);   init_v_or_w(m.w, 157,  13.0,   0.0);
        init_v_or_w(m.v,  16, -40.0,  40.0);   init_v_or_w(m.w,  16,  -1.0,   0.0);
        init_v_or_w(m.v, 158,  40.0, -40.0);   init_v_or_w(m.w, 158,  -1.0,   0.0);
        init_v_or_w(m.v,  15,   0.0,   0.0);
        init_v_or_w(m.w, 160,   0.0,   0.0);

        // 9. + 10. Polar cells — ventilate the polar caps (j=0..14 / 166..180), which were
        // previously left with ZERO meridional circulation (stagnant), a suspected driver of
        // the polar temperature/pressure build-up and the long-run pole blow-up. Two thin
        // overturning cells per cap continue the banded pattern to the pole (mirrors ATOM's
        // polar-cell construction), with the branch signs alternating from the 8. Ferrel cell
        // and zero meridional flow enforced at the poles j=0 / j=180. Magnitudes follow the
        // existing +-40 m/s cell pattern; tune as the physics dictates.
        // ONE polar cell per cap (rising centre), mirror-positioned N<->S. Unit:
        // [Z boundary] R(v:trop+40,surf-40 equatorward flank) P(rising centre, v=0)
        // R(v:trop-40,surf+40 poleward flank) [Z pole] — flow converges at the surface,
        // rises at the centre, diverges aloft (one closed rotation, no redundant divider).
        // northern cap: single rising cell — centre j=8, flanks j=12/4, boundaries j=15/0.
        init_v_or_w(m.v,  12,  40.0, -40.0);   init_v_or_w(m.w,  12,   0.0,   0.0);   // equatorward flank
        init_v_or_w(m.v,   8,   0.0,   0.0);   init_v_or_w(m.w,   8,   0.0,   0.0);   // rising centre
        init_v_or_w(m.v,   4, -40.0,  40.0);   init_v_or_w(m.w,   4,   0.0,   0.0);   // poleward flank
        init_v_or_w(m.v,   0,   0.0,   0.0);   init_v_or_w(m.w,   0,   0.0,   0.0);   // north pole
        // southern cap: kept as two cells (Jupiter is not equator-symmetric; not mirrored to N)
        init_v_or_w(m.v, 166, -40.0,  40.0);   init_v_or_w(m.w, 166,   0.0,   0.0);   //  9. equatorward flank
        init_v_or_w(m.v, 168,   0.0,   0.0);   init_v_or_w(m.w, 168,   0.0,   0.0);   //  9. rising centre
        init_v_or_w(m.v, 170,  40.0, -40.0);   init_v_or_w(m.w, 170,   0.0,   0.0);   //  9. poleward flank
        init_v_or_w(m.v, 172,   0.0,   0.0);   init_v_or_w(m.w, 172,   0.0,   0.0);   //  9./10. boundary
        init_v_or_w(m.v, 174, -40.0,  40.0);   init_v_or_w(m.w, 174,   0.0,   0.0);   // 10. equatorward flank
        init_v_or_w(m.v, 176,   0.0,   0.0);   init_v_or_w(m.w, 176,   0.0,   0.0);   // 10. rising centre
        init_v_or_w(m.v, 178,  40.0, -40.0);   init_v_or_w(m.w, 178,   0.0,   0.0);   // 10. poleward flank
        init_v_or_w(m.v, 180,   0.0,   0.0);   init_v_or_w(m.w, 180,   0.0,   0.0);   // south pole

        // --------------------------------------------------------------------
        // Interpolate between seeded latitudes — northern hemisphere
        // --------------------------------------------------------------------
        // u
        form_diagonals(m.u,   0,   4);  form_diagonals(m.u,   4,   8);   // polar cell
        form_diagonals(m.u,   8,  12);  form_diagonals(m.u,  12,  15);   // polar cell
        form_diagonals(m.u,  15,  16);  form_diagonals(m.u,  16,  17);
        form_diagonals(m.u,  17,  19);  form_diagonals(m.u,  19,  29);   // merged downwelling

        // (merged into 21->27)
        form_diagonals(m.u,  29,  30);
        form_diagonals(m.u,  30,  34);  form_diagonals(m.u,  34,  35);
        form_diagonals(m.u,  35,  37);  form_diagonals(m.u,  37,  40);
        form_diagonals(m.u,  40,  43);  form_diagonals(m.u,  43,  46);
        form_diagonals(m.u,  46,  48);  form_diagonals(m.u,  48,  49);
        form_diagonals(m.u,  49,  51);  form_diagonals(m.u,  51,  52);
        form_diagonals(m.u,  52,  53);  form_diagonals(m.u,  53,  54);
        form_diagonals(m.u,  54,  57);  form_diagonals(m.u,  57,  60);
        form_diagonals(m.u,  60,  64);  form_diagonals(m.u,  64,  66);
        form_diagonals(m.u,  66,  69);  form_diagonals(m.u,  69,  76);
        form_diagonals(m.u,  76,  83);  form_diagonals(m.u,  83,  90);
        // v
        form_diagonals(m.v,   0,   4);  form_diagonals(m.v,   4,   8);   // polar cell
        form_diagonals(m.v,   8,  12);  form_diagonals(m.v,  12,  15);   // polar cell
        form_diagonals(m.v,  15,  16);  form_diagonals(m.v,  16,  17);
        form_diagonals(m.v,  17,  19);  form_diagonals(m.v,  19,  29);   // merged downwelling

        // (merged into 21->27)
        form_diagonals(m.v,  29,  30);
        form_diagonals(m.v,  30,  34);  form_diagonals(m.v,  34,  35);
        form_diagonals(m.v,  35,  37);  form_diagonals(m.v,  37,  40);
        form_diagonals(m.v,  40,  43);  form_diagonals(m.v,  43,  46);
        form_diagonals(m.v,  46,  48);  form_diagonals(m.v,  48,  49);
        form_diagonals(m.v,  49,  51);  form_diagonals(m.v,  51,  52);
        form_diagonals(m.v,  52,  53);  form_diagonals(m.v,  53,  54);
        form_diagonals(m.v,  54,  57);  form_diagonals(m.v,  57,  60);
        form_diagonals(m.v,  60,  64);  form_diagonals(m.v,  64,  66);
        form_diagonals(m.v,  66,  69);  form_diagonals(m.v,  69,  76);
        form_diagonals(m.v,  76,  83);  form_diagonals(m.v,  83,  90);
        // w
        form_diagonals(m.w,   0,   4);  form_diagonals(m.w,   4,   8);   // polar cell
        form_diagonals(m.w,   8,  12);  form_diagonals(m.w,  12,  15);   // polar cell
        form_diagonals(m.w,  15,  16);  form_diagonals(m.w,  16,  17);
        form_diagonals(m.w,  17,  19);  form_diagonals(m.w,  19,  29);   // merged downwelling

        // (merged into 21->27)
        form_diagonals(m.w,  29,  30);
        form_diagonals(m.w,  30,  34);  form_diagonals(m.w,  34,  35);
        form_diagonals(m.w,  35,  37);  form_diagonals(m.w,  37,  40);
        form_diagonals(m.w,  40,  43);  form_diagonals(m.w,  43,  46);
        form_diagonals(m.w,  46,  48);  form_diagonals(m.w,  48,  49);
        form_diagonals(m.w,  49,  51);  form_diagonals(m.w,  51,  52);
        form_diagonals(m.w,  52,  53);  form_diagonals(m.w,  53,  54);
        form_diagonals(m.w,  54,  57);  form_diagonals(m.w,  57,  60);
        form_diagonals(m.w,  60,  64);  form_diagonals(m.w,  64,  66);
        form_diagonals(m.w,  66,  69);  form_diagonals(m.w,  69,  76);
        form_diagonals(m.w,  76,  83);  form_diagonals(m.w,  83,  90);

        // --------------------------------------------------------------------
        // Interpolate between seeded latitudes — southern hemisphere
        // --------------------------------------------------------------------
        // u
        form_diagonals(m.u,  90,  97);  form_diagonals(m.u,  97, 102);
        form_diagonals(m.u, 102, 107);  form_diagonals(m.u, 107, 111);
        form_diagonals(m.u, 111, 115);  form_diagonals(m.u, 115, 118);
        form_diagonals(m.u, 118, 121);  form_diagonals(m.u, 121, 123);
        form_diagonals(m.u, 123, 124);  form_diagonals(m.u, 124, 125);
        form_diagonals(m.u, 125, 126);  form_diagonals(m.u, 126, 127);
        form_diagonals(m.u, 127, 128);  form_diagonals(m.u, 128, 129);
        form_diagonals(m.u, 129, 130);  form_diagonals(m.u, 130, 136);
        form_diagonals(m.u, 136, 138);  form_diagonals(m.u, 138, 140);
        form_diagonals(m.u, 140, 142);  form_diagonals(m.u, 142, 143);
        form_diagonals(m.u, 143, 145);  form_diagonals(m.u, 145, 146);
        form_diagonals(m.u, 146, 149);  form_diagonals(m.u, 149, 152);
        form_diagonals(m.u, 152, 153);  form_diagonals(m.u, 153, 154);
        form_diagonals(m.u, 154, 155);  form_diagonals(m.u, 155, 157);
        form_diagonals(m.u, 157, 158);  form_diagonals(m.u, 158, 160);
        form_diagonals(m.u, 160, 166);  form_diagonals(m.u, 166, 168);   // polar cells
        form_diagonals(m.u, 168, 170);  form_diagonals(m.u, 170, 172);   // polar cells
        form_diagonals(m.u, 172, 174);  form_diagonals(m.u, 174, 176);   // polar cells
        form_diagonals(m.u, 176, 178);  form_diagonals(m.u, 178, 180);   // polar cells
        // v
        form_diagonals(m.v,  90,  97);  form_diagonals(m.v,  97, 104);
        form_diagonals(m.v, 102, 107);  form_diagonals(m.v, 107, 111);
        form_diagonals(m.v, 111, 115);  form_diagonals(m.v, 115, 118);
        form_diagonals(m.v, 118, 121);  form_diagonals(m.v, 121, 123);
        form_diagonals(m.v, 123, 124);  form_diagonals(m.v, 124, 125);
        form_diagonals(m.v, 125, 126);  form_diagonals(m.v, 126, 127);
        form_diagonals(m.v, 127, 128);  form_diagonals(m.v, 128, 129);
        form_diagonals(m.v, 129, 130);  form_diagonals(m.v, 130, 136);
        form_diagonals(m.v, 136, 138);  form_diagonals(m.v, 138, 140);
        form_diagonals(m.v, 140, 142);  form_diagonals(m.v, 142, 143);
        form_diagonals(m.v, 143, 145);  form_diagonals(m.v, 145, 146);
        form_diagonals(m.v, 146, 149);  form_diagonals(m.v, 149, 152);
        form_diagonals(m.v, 152, 153);  form_diagonals(m.v, 153, 154);
        form_diagonals(m.v, 154, 155);  form_diagonals(m.v, 155, 157);
        form_diagonals(m.v, 157, 158);  form_diagonals(m.v, 158, 160);
        form_diagonals(m.v, 160, 166);  form_diagonals(m.v, 166, 168);   // polar cells
        form_diagonals(m.v, 168, 170);  form_diagonals(m.v, 170, 172);   // polar cells
        form_diagonals(m.v, 172, 174);  form_diagonals(m.v, 174, 176);   // polar cells
        form_diagonals(m.v, 176, 178);  form_diagonals(m.v, 178, 180);   // polar cells
        // w
        form_diagonals(m.w,  90,  97);  form_diagonals(m.w,  97, 104);
        form_diagonals(m.w, 102, 107);  form_diagonals(m.w, 107, 111);
        form_diagonals(m.w, 111, 115);  form_diagonals(m.w, 115, 118);
        form_diagonals(m.w, 118, 121);  form_diagonals(m.w, 121, 123);
        form_diagonals(m.w, 123, 124);  form_diagonals(m.w, 124, 125);
        form_diagonals(m.w, 125, 126);  form_diagonals(m.w, 126, 127);
        form_diagonals(m.w, 127, 128);  form_diagonals(m.w, 128, 129);
        form_diagonals(m.w, 129, 130);  form_diagonals(m.w, 130, 136);
        form_diagonals(m.w, 136, 138);  form_diagonals(m.w, 138, 140);
        form_diagonals(m.w, 140, 142);  form_diagonals(m.w, 142, 143);
        form_diagonals(m.w, 143, 145);  form_diagonals(m.w, 145, 146);
        form_diagonals(m.w, 146, 149);  form_diagonals(m.w, 149, 152);
        form_diagonals(m.w, 152, 153);  form_diagonals(m.w, 153, 154);
        form_diagonals(m.w, 154, 155);  form_diagonals(m.w, 155, 157);
        form_diagonals(m.w, 157, 158);  form_diagonals(m.w, 158, 160);
        form_diagonals(m.w, 160, 166);  form_diagonals(m.w, 166, 168);   // polar cells (w overwritten by N->S mirror below)
        form_diagonals(m.w, 168, 170);  form_diagonals(m.w, 170, 172);   // polar cells
        form_diagonals(m.w, 172, 174);  form_diagonals(m.w, 174, 176);   // polar cells
        form_diagonals(m.w, 176, 178);  form_diagonals(m.w, 178, 180);   // polar cells

        // --------------------------------------------------------------------
        // Mirror w field north-south (Jupiter's hemispherical symmetry)
        // --------------------------------------------------------------------
        const int jm = m.jm;
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int k = 0; k < m.km; k++) {
                std::vector<double> tmp(jm);
                for (int j = 0; j < jm; j++) tmp[j] = m.w.x[i][j][k];
                for (int j = 0; j < jm; j++) m.w.x[i][j][k] = tmp[jm - 1 - j];
            }
        }

        // --------------------------------------------------------------------
        // Zero terrain cells; non-dimensionalise fluid cells — fused pass
        // --------------------------------------------------------------------
        const double inv_u_0 = 1.0 / m.u_0;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int k = 0; k < m.km; k++) {
                for (int j = 0; j < m.jm; j++) {
                    if (m.SeaMount.x[i][j][k] == 1.0) {
                        m.u.x[i][j][k] = 0.0;
                        m.v.x[i][j][k] = 0.0;
                        m.w.x[i][j][k] = 0.0;
                    } else {
                        m.u.x[i][j][k] *= inv_u_0;
                        m.v.x[i][j][k] *= inv_u_0;
                        m.w.x[i][j][k] *= inv_u_0;
                    }
                }
            }
        }

        // --------------------------------------------------------------------
        // The six latitude bands the integrator never reaches (ATJUP_ZERO_GHOST_BANDS)
        // --------------------------------------------------------------------
        // RungeKuttaJup runs j = 3 .. jm-4, so j = 0,1,2 and j = jm-3..jm-1 are never advanced:
        // whatever is written here stays, bit for bit, to the end of the run. Measured at
        // iteration 450 of a 6.4-hour run: j = 0,1,2 held -10.0, -10.0, -20.0 m/s of RADIAL
        // velocity and j = 178,179,180 held -40.0, -20.0, -20.0, all identical to iteration 0,
        // while the largest |u| anywhere in the integrated range was 7.1 m/s. They are therefore
        // not initial data at all but a permanent boundary condition on the computed domain, five
        // times stronger than the solution it borders, and nothing in the model can relax them.
        // ChemMassRateJup already zeroes the same bands for its own reasons: without a dynamical
        // sink the quadratic reaction term there grows without bound.
        //
        // Setting them to zero is a diagnostic, not a repair — the honest repair is either to
        // integrate them or to give them a proper polar boundary condition. This knob measures
        // what that prescribed forcing is worth. Default 0 = every existing run bit-identical.
        static const int zero_bands = [](){
            const char* e = getenv("ATJUP_ZERO_GHOST_BANDS"); return e ? atoi(e) : 0; }();
        if(zero_bands){
            printf("      ATJUP: ATJUP_ZERO_GHOST_BANDS - zeroing u,v,w in j=0..2 and j=%d..%d\n",
                   m.jm - 3, m.jm - 1);
            #pragma omp parallel for collapse(2) schedule(static)
            for (int i = 0; i < m.im; i++) {
                for (int k = 0; k < m.km; k++) {
                    for (int b = 0; b < 6; b++) {
                        const int j = (b < 3) ? b : m.jm - 6 + b;
                        m.u.x[i][j][k] = 0.0;
                        m.v.x[i][j][k] = 0.0;
                        m.w.x[i][j][k] = 0.0;
                    }
                }
            }
        }

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" Time measured: %.3f seconds for VelocityInitializerJup\n",
               elapsed.count() * 1e-9);
        std::cout << "      ATJUP: VelocityInitializerJup::compute ended" << std::endl;
    }

private:
    cJupiterModel& m;

    // ========================================================================
    // Amplitude of the RADIAL branch of the overturning cells at latitude j.
    //
    // The names below (P/R, "prograde"/"retrograde") and the magnitude 80 m/s are a zonal-jet
    // vocabulary, but init_u() writes them into u, which is the radial component
    // (cJupiterModel.h:898). The zonal jet is w and is seeded separately in compute(), with
    // 13..162 m/s. The ancestor this construction came from, ATOM's
    // VelocityInitializer::init_u, uses 0.0116..0.0289 m/s here — three orders of magnitude
    // less — so the 80 is very likely a zonal number on the wrong component.
    //
    // What the grid allows: a cell is 1216 km zonally against 3.5 km radially, so continuity
    // (du/dr = -div_h v_h) caps the radial velocity at ~11 m/s for a structure one cell wide
    // and at ~1 m/s for a cell spanning several degrees of latitude. The realised peak here is
    // ratio 0.5 x 80 = 40 m/s, above both.
    //
    // ATJUP_U_INIT_SCALE multiplies the amplitude so the two can be separated by measurement
    // without committing to either reading. Default 1.0 = unchanged.
    // ========================================================================
    static double u_init_scale()
    {
        static const double v = [](){
            const char* e = getenv("ATJUP_U_INIT_SCALE"); return e ? atof(e) : 1.0; }();
        return v;
    }

    static double u_amplitude(int j)
    {
        const double s = u_init_scale();
        const double P =  80.0 * s;   // prograde jet
        const double R = -80.0 * s;   // retrograde jet
        const double Z =   0.0;       // cell centre (node)

        switch (j) {
            // northern hemisphere — 8 cell pairs
            case  15: return Z;  case  16: return R;
            case  17: return P;  // 8. Ferrel/Hadley centre
            case  19: return R;  // j21..27 boundary removed (gap closed)
            // case 24 (P centre) deleted — redundant updraft; interpolated as sink
            case  29: return R;
            case  30: return P;  // 6. Ferrel/Hadley centre
            case  34: return R;  case  35: return Z;  case  37: return R;
            case  40: return P;  // 5. Ferrel/Hadley centre
            case  43: return R;  case  46: return Z;  case  48: return R;
            case  49: return P;  // 4. Ferrel/Hadley centre
            case  51: return R;  case  52: return Z;  case  53: return R;
            case  54: return P;  // 3. Ferrel/Hadley centre
            case  57: return R;  case  60: return Z;  case  64: return R;
            case  66: return P;  // 2. Ferrel/Hadley centre
            case  69: return R;  case  76: return Z;  case  83: return R;
            // equator
            case  90: return P;
            // southern hemisphere — mirror of north
            case  97: return R;  case 102: return Z;  case 107: return R;
            case 111: return P;  // 1. Hadley/Ferrel centre
            case 115: return R;  case 118: return Z;  case 121: return R;
            case 123: return P;  // 2. Hadley/Ferrel centre
            case 124: return R;  case 125: return Z;  case 126: return R;
            case 127: return P;  // 3. Hadley/Ferrel centre
            case 128: return R;  case 129: return Z;  case 130: return R;
            case 136: return P;  // 4. Hadley/Ferrel centre
            case 138: return R;  case 140: return Z;  case 142: return R;
            case 143: return P;  // 5. Hadley/Ferrel centre
            case 145: return R;  case 146: return Z;  case 149: return R;
            case 152: return P;  // 6. Hadley/Ferrel centre
            case 153: return R;  case 154: return Z;  case 155: return R;
            case 157: return P;  // 7. Hadley/Ferrel centre
            case 158: return R;  case 160: return Z;
            // 9./10. polar cells (cap ventilation) — give the polar overturning a radial
            // (vertical) branch, so the cells are not flat. Boundaries/poles = Z (node),
            // branches alternate P/R like the neighbouring cells.
            case  12: return R;  case   8: return P;  case   4: return R;  case   0: return Z;  // north polar cell (centre j8)
            case 166: return R;  case 168: return P;  case 170: return R;  case 172: return Z;  // south polar cell 9
            case 174: return R;  case 176: return P;  case 178: return R;  case 180: return Z;  // south polar cell 10 + pole
            default:  return Z;
        }
    }

    // ========================================================================
    // Vertical u-profile: triangular ramp up to tropopause
    //   ascending  (h <  half_h): ratio = h / (half_h / 0.5)
    //   descending (h >= half_h): ratio = (tropo_h - h) / half_h
    // ========================================================================
    void init_u(Array& u, int j)
    {
        const double amp            = u_amplitude(j);
        const int    tl             = m.im_tropopause[j];
        const double tropo_h        = m.get_layer_height(tl);
        const double half_h         = tropo_h / 1.5;
        const double inv_ascent     = 0.5 / half_h;   // = 1 / (half_h/0.5)
        const double inv_descent    = 1.0 / half_h;

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < tl; i++) {
                const double h     = m.get_layer_height(i);
                const double ratio = (h < half_h)
                    ? h * inv_ascent
                    : (tropo_h - h) * inv_descent;
                u.x[i][j][k] = amp * ratio;
            }
        }
    }

    // ========================================================================
    // Linear v/w profile from surface value to tropopause value, then decay
    // ========================================================================
    void init_v_or_w(Array& v_or_w, int j, double coeff_trop, double coeff_sl)
    {
        const int    tl          = m.im_tropopause[j];
        const double tropo_h     = m.get_layer_height(tl);
        const double inv_tropo_h = 1.0 / tropo_h;

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            double sl = coeff_sl;
            if (is_ocean_surface(m.SeaMount, 20, j, k))
                sl = v_or_w.x[0][j][k];

            const double slope = (coeff_trop - sl) * inv_tropo_h;
            for (int i = 0; i < tl; i++) {
                v_or_w.x[i][j][k] = slope * m.get_layer_height(i) + sl;
            }
        }

        init_v_or_w_above_tropopause(v_or_w, j, coeff_trop);
    }

    void init_v_or_w_above_tropopause(Array& v_or_w, int j, double coeff)
    {
        const int tl = m.im_tropopause[j];
        if (tl >= m.im - 1) return;

        const double h_top     = m.get_layer_height(m.im - 1);
        const double inv_range = 1.0 / (h_top - m.get_layer_height(tl));

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = tl; i < m.im; i++) {
                v_or_w.x[i][j][k] = coeff * (h_top - m.get_layer_height(i)) * inv_range;
            }
        }
    }

    // ========================================================================
    // Linear interpolation across j in [start, end) — one i-level at a time
    // ========================================================================
    void form_diagonals(Array& a, int start, int end)
    {
        const double inv_range = 1.0 / (double)(end - start);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                const double a_start = a.x[i][start][k];
                const double slope   = (a.x[i][end][k] - a_start) * inv_range;
                for (int j = start; j < end; j++) {
                    a.x[i][j][k] = slope * (double)(j - start) + a_start;
                }
            }
        }
    }
};
