# ATJUP — notes for Claude Code

## Optional modules: knob tables

Moved verbatim from README.md on 2026-11-26 (the README keeps a short summary of this section). A new knob gets its row HERE, not in the README.

The three modules run by default, and so do their **coupling** knobs: `ATJUP_RAD_COUPLING` 1
(since 2026-10-27), `ATJUP_PRECIP_COUPLING` 1 (since 2026-07-31) and `ATJUP_TURB_COUPLING` 1 (since
2026-11-01), so a stock run feeds `Q_rad`, `Q_precip` and the eddy viscosity back into the
temperature, momentum and species equations. Set a module's variable to 0 to switch it off, or its
coupling knob to 0 to keep it diagnostic-only.

**Radiation**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_RADIATION` | 1 | grey two-stream solve; fills `Q_rad`, `radiation`, `epsilon` |
| `ATJUP_RAD_COUPLING` | 1 | feed `Q_rad` into the temperature equation (see note below); default 1 since 2026-10-27, 0 = off |
| `ATJUP_SOLAR` | 1 | absorbed shortwave channel (only acts with `ATJUP_RADIATION`) |
| `ATJUP_CIA_STRENGTH` | 3.9 | scale the H₂/He collision-induced opacity; default via `rad_default_cia_strength()` since 2026-10-27 (photosphere ~0.26 bar, OLR 13.7 W/m² on the corrected composition), 1 = the old value |
| `ATJUP_RAD_BOTTOM` | 0 | 1 = no radiative divergence in the bottom layer (`Fu = 2 B - Fd` above it; removes the level-0 jump in Radiation / Q_rad). Opt-in on ATJUP: the bottom is not opaque there, so it changes the state (with RAD_COUPLING 1); default via `rad_default_bottom()` (row GIANTS-RAD-BOTTOM, 2026-11-01) |
| `ATJUP_OPACITY_STRENGTH` | 1 | scale the gas-band and cloud opacity. The shared `opac_cal` 0.025 was tuned on the old fields (CH4 ~18x and the NH3 ice deck ~1000x the corrected state); on the corrected state the band + cloud opacity barely acts (OPACITY_STRENGTH ~650 would be needed for 0.3 bar, CP1_1000 sweep 10-26). Decided 2026-10-27: it stays at 1, `opac_cal` stays shared at 0.025, and `ATJUP_CIA_STRENGTH` 3.9 carries the photosphere -- a large factor here would tie the radiation to the still-forming NH3 deck |

The shared `Radiation.h` also prints, every run with `ATJUP_RADIATION=1`, the τ=1 photosphere
pressure **and the temperature there** beside the blackbody flux it implies. That second number is
what separates "the scheme emits wrongly" from "the scheme was handed a column that is too warm" —
`T(τ=1) − T_eff(OLR)` is the scheme's own excess, `T(τ=1) − T_eff(in)` is the column's. ATJUP is the
only one of the four whose photosphere temperature is close to right, and the only one that has
never needed a viscosity or pressure constant corrected. `epsilon` reaches ParaView as `Emissivity`
alongside `Radiation` and `Q_rad_mW_m3`, in all four views.

ATJUP carries **no `thermalmassflux` term in `rhs_t`** — it has `radiation_t + precip_t` in that slot
instead — so it has no `ATJUP_THERMAL_MASSFLUX` knob. The other three models do.

**Microphysics**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_PRECIP` | 1 | three-category microphysics + NH₄SH sedimentation |
| `ATJUP_PRECIP_COUPLING` | 1 (since 2026-07-31) | feed the latent heat `Q_precip` into the temperature equation, together with the condensate depletion it reports (`S_precip_*`); 0 = diagnostic only |
| `ATJUP_CH4_LOCAL_RHO` | 1 (since 2026-10-26) | CH4 as a fixed mass fraction of the local density, at init and after a restart load; 0 restores the constant-density field (0.19 kg/m³ everywhere, a mass fraction above 1 over ~80 km) |
| `ATJUP_T_DIFF_SCALE` | 0.01 (since 2026-10-26) | factor on the laminar heat diffusivity 1/(re pr) of the temperature equation (the eddy part is untouched), the counterpart of `ATJUP_SPECIES_DIFF_SCALE`; at 1 (the old value, ~2e4 m²/s) the corrected state drifts -11 K at the bottom and +5.5 K aloft per 7000 s and the NH3 ice deck evaporates |
| `ATJUP_COMPOSITION` | 1 (since 2026-10-26) | the observed composition: mass shares from He/H₂ 0.157, CH₄/H₂ 2.37·10⁻³, H₂S/H₂ 8.9·10⁻⁵ (Galileo), H₂O at 2.5× solar O/H, NH₃ 362 ppm (Juno), r_mix held; ep from the mixture, p_stat with R_mix; H₂O / NH₃ / H₂S start as ρ_local × min(deep share, cold trap), no condensate. From scratch only — on a restart H₂O, H₂S and p_stat keep the file's state |
| `ATJUP_INIT_LAPSE_ADIABAT=<f>` | 0 | initial tropospheric lapse rate gam = f x g/cp_mix (2.1946 K/km at f = 1) instead of the config's 2.0 K/km, which starts the deep atmosphere 0.2 K/km subadiabatic (N² ~2e-5). The profile pivots about the 1-bar level (T = t_ref, where p_stat = p_ref), so the 1-bar height and p from about i20 up stay within 1 %; the deep column warms (+15.6 K at the bottom, 12.73 -> 11.91 bar at the equator). Overrides the member gam, so the p_stat exponent and SaturationAdjustment follow; start such arms from scratch. Opt-in (2026-10-29; 2-it test A2 vs A0: N² ~0 in i3-28; no dynamics arm, LA0 stopped at it 260) |
| `ATJUP_NH4SH_EQ` | 1 (since 2026-10-26) | NH₄SH as an equilibrium adjustment on the Lewis K_p (log₁₀ K_p = 14.82 − 4705/T, atm²) with its heat of formation (1.76·10⁶ J/kg), instead of the rate law; parcel-swept initial cloud |
| `ATJUP_SVP_FIT` | 1 (since 2026-10-26) | NH₃ liquid / NH₃ ice / H₂O ice vapour pressures from NIST Antoine, Fray & Schmitt (2009) and Murphy & Koop (2005), fitted into the model's four-coefficient form (coefficients in `setSvpFit`) |

**Turbulence**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_TURB` | *unset* | override the switch that `turb_model` sets: 0 = closure off, anything else = on |
| `ATJUP_TURB_MODEL` | *param* | override `turb_model` (`none`, `k_epsilon`, `k_omega`, `k_omega_SST`) |
| `ATJUP_TURB_COUPLING` | 1 | feed the eddy viscosity into momentum, heat and species diffusion; default 1 since 2026-11-01, 0 restores |
| `ATJUP_TURB_FREE` | 1 | no wall law and no surface boundary layer in the closure: uniform background k / omega at the start, u_tau the calm-air floor everywhere, the k-epsilon / SST wall distance 'far' (free-flow form); 0 restores the ATOM surface picture (row GIANTS-TURB-FREE, 2026-11-01). ATJUP_TURB_WALL_HOLD / _WALL_OMEGA act only above an obstacle, so with NO_SEAMOUNT 1 they do nothing |
| `ATJUP_TURB_CLAMP_PASS` | 1 | the closure's clamp (k >= 0, omega >= min) runs as a pass before its source loop and the k floor / caps as a pass after it, so the loop no longer reads half-updated neighbours and the result does not depend on the thread count; 0 restores the in-loop form (row SAT-TURB-THREADS, 2026-11-01) |
| `ATJUP_TURB_BUOY` | 0 (off) | buoyancy term of the k equation, -nu_t N^2 / Pr_t (Pr_t 0.9), which the shear-only closure lacked; 1 = dry N^2 = (g/T)(dT/dz + g/cp), stable layers only; 2 = also N^2 < 0 as production; 3 = 1 with the moist N^2 (Durran & Klemp 1982) in cells holding condensate, gases from the model's turb_moist_gas() overload (ATSAT only so far; without one =3 is =1). Omega untouched (C3 = 0). Measured on ATSAT: the eddy heat flux nue_t/Pr_t had warmed the top +0.7 K / 100 it; =1 and =3 stop it, but k then decays to the k_bg floor above the lowest levels (=3 ~800 iterations later than =1) (row SAT-LID-CUBIC, 2026-11-06) |
| `ATJUP_PRESS_BAROTROPIC` | 0 (off) | column-mean (barotropic) mode of p_dyn solved EXACTLY after every sweep of the shared solver: a column-uniform q(theta, phi) zeroes the column-summed residual (DFT in phi, tridiagonal in theta per mode, theta ends folded as the boundary pass sets them). WHY: on a thin shell the radial stencil weight swamps the horizontal (ATSAT 1600 vs 0.24), so one pointwise sweep per call never moves p_dyn's horizontal column mean (~2e8 sweeps pole to pole); PRESS_LINE does not reach it either. Needs BC_PHI_PERIODIC 1 and no solid body (else skipped with one notice). Measured on ATSAT (B3L 4500 -> 5000, re 1000, satchk/pressbaro_2026-11-07 P1 vs retest_2026-11-06 R3): converges in 3 calls; the solved meridional gradient is within 2-4 % of the frozen initGeoBalance one, v max 115.8 -> 113.4 m/s, OLR / residuum unchanged -- the frozen p_dyn was NOT what drives the meridional cells (row SAT-SPINUP-DRIFT, 2026-11-07) |
| `ATJUP_NUE_MAX` | 1e5 | eddy-viscosity ceiling [m²/s] — a runaway guard, not the operative limiter |
| `ATJUP_TURB_BUOY` | 0 | buoyancy production in the k equation, G_b = -(nu_t/Pr_t) N², N² = (g/T)(dT/dz + g/cp_mix) (dry); the ω equation takes it on the unstable side only (C3 = 1 / 0). Pr_t = `ATJUP_TURB_BUOY_PRT` (0.9, ATOM's pr_turb). Opt-in (row JUP-FREE-GRS, 2026-10-29): on the REF / LR6 state the whole troposphere is stable (N² 2e-5 deep to 5e-4 aloft, Ri 50-900), so the offline estimate on LR6 5000 has k decaying in ~100 % of cells; no dynamics arm run (LA1 stopped at it 260) |
| `ATJUP_ABL_TAPER` | 0 | restore ATOM's boundary-layer taper of the eddy viscosity |
| `ATJUP_TURB_WALL_HOLD` | 1 | reassert the closure's wall condition above the obstacle after every Runge-Kutta step, as ATOM does; with 0 (the behaviour until 2026-10-09) the step overwrites it and k there runs to its ceiling (1000 m²/s² at iteration 368) |
| `ATJUP_TURB_WALL_OMEGA` | 1 | give the wall cell above the obstacle a log-layer wall function for ω; 0 (the behaviour until 2026-10-09) leaves the viscous-sublayer value, which on a 3.5 km layer is zero and drains ω from the shear layer above it. Max k at iteration 896: 502 with 0, 231 m²/s² with 1, levelling off |

Worth stating plainly, because it is easy to believe otherwise: **`turb_model` in the XML does not
Since 2026-10-09 **`turb_model` is the switch of the closure**, as in ATOM: `none` (or `laminar`,
or empty) is off, any of the three models is on, and "on" means both halves — `TurbulenceJup`
fills `tke`, `dis` and `nue`, and k and ω are advanced by the Runge-Kutta step. The parameter
defaults to `k_omega_SST`, so a stock run has the closure on. `ATJUP_TURB=0` still switches it
off from the environment, whatever the XML says. Without `ATJUP_TURB_COUPLING=1` the eddy
viscosity is computed but never reaches the momentum, heat or species equations; that is the one
place ATJUP still differs from ATOM, where it always does.

Two earlier states, for reading old runs. Until `8df7626` `ATJUP_TURB` defaulted to 0 and
`turb_model` only selected the model, so a run "with k-ω SST" that set only the XML had no
turbulence in it. From `8df7626` until 2026-10-09 the variable was read in two places with
different defaults (1 in `cJupiterModel.cpp`, 0 in `RHS_Jup_Turb.cpp`), so a run that left it
unset computed the eddy viscosity but did not advance k and ω: `tke` stayed at its initial
value. `ATJUP_TURB=1` set explicitly gave then what the default gives now, byte for byte.

**Initial and boundary conditions**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_IC_TROPO_CLAMP` | 1 | tropopause clamp on the temperature initial condition |
| `ATJUP_T_TROPOPAUSE_MIN` | 110 | tropopause minimum temperature [K] |
| `ATJUP_VEL_SHAPIRO_ORDER` | 2 | 2 = 1-2-1 filter, 4 = shear-preserving |
| `ATJUP_BC_RIGID_LID` | 1 | u = 0 at both radial walls |
| `ATJUP_BC_PHI_PERIODIC` | 1 | since 2026-10-13 (0 restores the reconstructed seam, byte-identical to the head before). The longitude grid carries the same meridian twice, at k = 0 and k = km-1, and no loop that differences in phi integrated it: `bcPhi` set it to the mean of two one-sided extrapolations, the pressure solver did the same to `p_dyn`, the Shapiro filter took the two copies for each other's neighbour, and the output routine reconstructed T at the seam in the state at every write. With 1 the meridian is an ordinary cell — Runge-Kutta, RHS, chemistry, latent heat, forces, pressure solver, turbulence closure, flux limiter and filter take their phi neighbours round the ring of the km-1 distinct meridians — and k = km-1 is its copy. Measured here at 8 iterations (satchk/giants, JPH0 against JPH8): 115 of 120 printed extrema identical, continuity residuum 0.035253 unchanged (it sits at the obstacle); nothing longer was run. `bcSolidGround` still works k = 1..km-2: the obstacle spans k = 140..220 and does not reach the seam |
| `ATJUP_METRIC_SPLIT` | 1 | since 2026-10-13 (0 restores, byte-identical to the head before): the radial metric term (v²+w²)/r of `rhs_u` is carried by the radial integral of itself, rebuilt every step, instead of by the pressure solver — ATOM_Precipitation's hydrostatic split applied to the other radial body force. `rhs_u` drops the term (it and the radial gradient of that pressure cancel by construction), `rhs_v` and `rhs_w` take its horizontal gradient, and `p_dyn` starts at zero. Why: the solver cannot hold the gradient that balances a radial body force (measured on ATNEPT with `ATNEPT_UB_DIAG`: after 224 iterations the balanced initial `p_dyn` had lost nine tenths of its gradient next to the radial walls; residuum there 0.008808 → 0.000341 with the split). The pressure is `p_hydro`, integrated over the fluid part of each column, as the buoyancy split does around the obstacle. Measured here at 8 iterations (satchk/giants, JPDd against JMDd): 100 of 120 printed extrema identical, max u 0.00133 → 0.00065, max `p_dyn` 0.00013 → 0.000012, continuity residuum 0.035253 → 0.035234 (it sits at the obstacle and is not this). 224 iterations, 0 against 1 (JM0_224, JM1_224): continuity residuum 0.030077 → 0.031376 (at the obstacle; from iteration ~112 it stops falling — not looked into), printed max u 0.0569 → 0.0142, max `p_dyn` 0.000078 → 0.000016, H2O column −0.0595 → −0.0596 % |
| `ATJUP_CORIOLIS_SPLIT` | 1 | since 2026-10-13 (0 restores, byte-identical to the head before; acts only with `ATJUP_METRIC_SPLIT`): the radial Coriolis term of `rhs_u` joins the metric term in the same column integral and leaves `rhs_u` — it is a radial body force too, and with the metric term gone it was what still drove u (ATNEPT, equatorial band: −0.014 m/s in 224 iterations at every level). It enters the integral with the true sine, not a metric floor. 224 iterations, 0 against 1 (satchk/giants/JM1_224, JCS224): continuity residuum 0.031376 → 0.031864 (it sits at the obstacle; 0.030077 without either split — not looked into), printed max u 0.0142 → 0.0156, H2O column −0.0596 → −0.0598 % |
| `ATJUP_PGRAD_WRAP` | 1 | since 2026-10-17 (0 restores, byte-identical to the head before): the φ neighbours of the pressure gradients (`pgrad` in `RHS_Jup_Turb.cpp`, for `p_dyn` and `p_hydro`) through `phi_wrap`, as every other φ derivative takes them. Since `ATJUP_BC_PHI_PERIODIC` integrates the seam meridian k = 0, `k − 1` read `P.x[i][j][-1]`, out of bounds — the last cell of row j−1, i.e. the MERIDIONAL difference of a zonally uniform field entering `rhs_w` as a zonal force on that one meridian. With the balanced p_dyn of `ATJUP_GEO_BALANCE` it drove w at k = 0 +8.6 m/s in 100 iterations at 86 N (satchk/geo_2026-10-17/GEO224_seambug); in the default run p_dyn is ~1e-5 and the fix changes w by at most 0.0004 m/s at k = 0 by iteration 200 and nothing else in the printed digits (PW224 against cor_2026-10-17/JCF1000). ATSAT / ATNEPT / ATURAN never had it |
| `ATJUP_GEO_BALANCE` | 1 (since 2026-10-22; was 0) | 0 restores the head before (f857a9c) byte for byte. Since 2026-10-17: the geostrophically balanced start of ATNEPT / ATURAN (from ATOM_Precipitation). Switches on, each with its own knob: `ATJUP_ND_COR` 1 (Coriolis × L/u₀ = 1400), `ATJUP_PDYN_GEO` 1, `ATJUP_PRESS_INCR` 1, `ATJUP_DAMP_PDYN` 0, `ATJUP_RHO_PDYN` 0. WHY: with `ATJUP_ND_COR=1` alone v grows at exactly f·w and nothing answers it (cor_2026-10-17, 1000 it: +2.66 m/s per 100 it at 25 N; dv/(w·it) = f·dt to 1 %). 224 iterations, 6 threads (satchk/geo_2026-10-17/GEO224b against cor_2026-10-17/JCF1000, JNC1000): largest v change against the default at 200, zonal slice 0.91 m/s (ND_COR alone 5.59), level 20 0.44 (3.36); residuum at 208 0.0331 (default 0.0319, ND_COR alone 0.0317), now at the obstacle 22 S 171 E 115 km; max u 0.085 m/s at the obstacle (default 0.015) — the balance is built on zonal means and the obstacle flank takes one-sided gradients. Left unbalanced: the prescribed meridional cells, which full Coriolis turns into w (+3.2 m/s near 86 N at the bottom by 200, as with ND_COR alone) |
| `ATJUP_CELL_HOLD` | follows `ATJUP_GEO_BALANCE` | since 2026-10-22: the zonal Coriolis component acts on the departure from the frozen initial zonal-mean meridional cells (v - vbar0, u - ubar0) -- a prescribed zonal force that holds the cells in rotating balance; (i, j) only, pins nothing in longitude. Without it full Coriolis turns the +-40 m/s cells into w (satchk/cor_2026-10-22 CR2 vs CR4) |
| `ATJUP_V_INIT_SCALE` | 1 | since 2026-10-22: multiplies the prescribed meridional cells; 0 starts v at zero (jets held, but the vertical motion and the ice-updraft link go with it, CR3) |
| `ATJUP_PDYN_GEO` | follows `ATJUP_GEO_BALANCE` | `initGeoBalance()` (InitValues_Jup.cpp), after the initial boundary passes: p_dyn = the barotropic pressure that balances Coriolis and metric force of each latitude's jet-level zonal-mean wind in `rhs_v`, integrated from the north pole by p[j+1] = p[j−1] + 2Δθ G[j] (the centred difference returns G exactly); `p_baro` (frozen, added to `p_hydro`) the per-level rest. Needs `ATJUP_METRIC_SPLIT`. p_dyn −1.25 … +1.02 bar, p_baro −16.3 … +1.5 (model unit) |
| `ATJUP_PRESS_INCR` | 0 (since 2026-10-25; followed `ATJUP_GEO_BALANCE` before, 1 restores) | WHY 0: with it the radial wind in the top fluid cell grew and flipped sign against the column = the 39 S lid residuum band (CR5 0.011 → 0.062 over 7000 s); satchk/pinc_2026-10-25 PI1 / PI0 (vortex + PRESS_LINE, 1000 it): residuum max 0.048 → 0.017, 39 S band 0.031 → 0.0005, vortex / jets / ice unchanged (|Δw| < 1 m/s). Before: the shared solver's defect-correction source (`planet/PressureSolver.h`, `<TAG>_PRESS_INCR`, default 0 there): Lc(p_new) = div(rhs) + Lc(p_old), so a balanced p_dyn is a fixed point; with a rigid lid the wall-normal source is zero at the radial walls |
| `ATJUP_DAMP_PDYN` | 0 with `ATJUP_GEO_BALANCE`, else 1 | the Shapiro pass on p_dyn after the solver; on ATNEPT it halved the balancing gradient in 224 iterations |
| `ATJUP_RHO_PDYN` | 0 with `ATJUP_GEO_BALANCE`, else 1 | `rho_mix` built with p_dyn (1) or from p_stat alone (0), as the buoyancy is. The balanced p_dyn is a potential on the constant r_mix; ±1 bar against 0.02 bar under the lid would make the density the chemistry reads negative |
| `ATJUP_UB_DIAG=<n>` | off | print-only (2026-10-13, `planet/UbJup.h`, ported from ATNEPT): every n iterations, which stage of the loop and which term of `rhs_u` changed the level-mean radial wind, for the whole level and for the equatorial band (`ATJUP_UB_BAND=<degrees>`, default 5), with the radial pressure force and what the pressure solver and the filter on `p_dyn` did to it. The level means include the cells inside the obstacle, the pressure force is the plain centred difference, and the filter on u runs only with `ATJUP_SHAPIRO_VEL_INLOOP`. Unset and set both leave every output file identical (8 iterations) |
| `ATJUP_UB_LAT=<deg N>` | 0 | print-only (2026-10-27): centre of the `ATJUP_UB_DIAG` band (0 = the equator as before), for the zonal lid bands of row JUP-RES23N |
| `ATJUP_BC_TOP_TAPER` | 3 | top sponge over the top three layers. 3 (default since 2026-10-29, row JUP-RES23N) = zonal-mean-preserving: v = vbar + f (v - vbar) with f = 2/3, 1/3 and the lid at the im-2 mean, so only the deviations from the zonal mean (vortex, waves) are damped; 1 = the old taper of v, w to zero; 0 = off. Measured with `ATJUP_W_UPPER_LID_FRAC=0.85` to 3000 (WS3, satchk/sponge_2026-10-29): i39 lid bands 0.001 (old taper 0.014-0.033), vortex 0 at i39, a local residuum peak 0.064 left at the vortex top |
| `ATJUP_PHYDRO_TOP_TAPER` | 0 | 1 = drop the horizontal p_hydro force in the two top fluid layers, where `ATJUP_BC_TOP_TAPER` holds v and w anyway; 2 = scale it 2/3, 1/3. Opt-in (row JUP-RES23N, 2026-10-27/28): over 1000 iterations from WK0 2000 it lowered the domain max of the lid residuum (0.075 -> 0.053) but raised the 53 S / 23 N bands by 78 / 46 % and the interior by 50 %, and moved the upper T by up to 11 K, so it stays off |
| `ATJUP_W_UPPER_LID_FRAC` | 0.85 | the initial zonal wind above the tropopause falls linearly to this fraction of its tropopause value at the lid (0 = to zero, the old start; 1 = no decay); v keeps the decay. Default since 2026-10-29 together with `ATJUP_BC_TOP_TAPER=3` (row JUP-RES23N, option 2, 2026-10-28): the decay to 0 needs about +-60 K of T contrast for thermal-wind balance, which the isothermal start does not have. Measured to 3000 (satchk/upper_2026-10-28): with the top taper on (WT1) it equals REF; with `ATJUP_BC_TOP_TAPER=0` (WT0, frac 0.85) the i39 lid bands drop 20-30x, but the jets stay barotropic to the lid and the vortex reaches it (20 cells on the 150 m/s cap) |
| `ATJUP_BC_POLE_COPY` | 1 | plain copy instead of extrapolation at the poles |
| `ATJUP_BC_POLE_FILL` | 1 | default since 2026-10-11, 0 restores. Also sets the two rows next to each pole (`j = 1, 2` and their mirrors) from the first row the Runge-Kutta step integrates (`j = 3`): scalars and `u` by copy, `v`, `w` and the flux fields by a linear ramp to their zero on the axis. Without it those rows are written by nobody and keep the initial condition for the whole run |
| `ATJUP_R_NH3_ADD` | *param* | deep well-mixed NH₃ density [kg/m³] |

**Geometry, scaling and the obstacle**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_METRIC_RADIUS` | 0 (off) | planetary radius [km] for the metric terms — see *Geometry* below |
| `ATJUP_SINTHE_MIN` | 0.55 | polar floor on sin θ; 0.55 means the floor is active poleward of 56.6° |
| `ATJUP_WALL_NUE` | 4 | wall eddy viscosity at the obstacle, in multiples of 1/re; 0 disables |
| `ATJUP_WALL_NUE_LAYERS` | 3 | its ramp depth in cells |
| `ATJUP_LOCAL_RHO` | 1 (since 2026-10-22) | use the local density instead of the constant r_mix in q_sat, latent heat and the Lewis groups |
| `ATJUP_NONDIM` | 0 | put the body forces into the units `rhs_u` is written in — see *Scaling* |
| `ATJUP_ND_BUOY` `_COR` `_CENT` | *follow `ATJUP_NONDIM`* | the three factors individually, for attribution |
| `ATJUP_ADIABATIC` | 1 (since 2026-10-22; was: follows `ATJUP_NONDIM`, i.e. off) | compression work against p_stat, i.e. the dry adiabatic lapse rate; 0 restores |
| `ATJUP_ADIABATIC_RHO` | 1 | that term divided by the local density instead of the constant r_mix (0 = old: -(g/cp) rho/r_mix, 0.22 instead of 2.07 K/km at the NH3 deck) |
| `ATJUP_MOUNT_T_AMBIENT` | 1 | every obstacle cell at the zonal-mean air temperature of its row (level, latitude): no temperature jump at the cone wall (0 = old: values extrapolated from the fluid side, cold sheath on the flanks) |
| `ATJUP_CONV_ADJ` | *follows `ATJUP_NONDIM`* | dry convective adjustment of superadiabatic columns |
| `ATJUP_CONV_ADJ_LAPSE` | 1 | scale its critical lapse rate (1 = the dry adiabat) |
| `ATJUP_CONV_ADJ_PASSES` | 64 | cap on sweeps per column; a warning is printed if it is reached |
| `ATJUP_MCA` | 1 (since 2026-10-26) | moist convective adjustment with the Ledoux criterion every N iterations, 0 = off (see below) |
| `ATJUP_MCA_LEDOUX` `_LOADING` `_RELAX` `_PASSES` | 1, 1, 1, 64 | composition in the density test (0 = temperature only), condensate weight, fraction applied per call, sweep cap |
| `ATJUP_MCA_DEBUG` | — | `j,k`: print the parcel test of one column on the first call |
| `ATJUP_MCA_INTERIOR` | 1 | 1 = leave the radial boundary planes i = 0 and i = im-1 out of the MCA column (shared `ConvectiveAdjustment.h`, `mca_default_interior()`; 1 by default on all four since 2026-10-30, on ATURAN since 10-28, where mixing the floor plane looped against `bcRadius`; 0 = the old column) |
| `ATJUP_MCA_SMOOTH=<n>` | 0 | stability tested on a column smoothed in latitude (<n> 1-2-1 passes over the rows j-1, j, j+1 of h_L, q_t and the tracers, snapshot before the column loop); the union of the segments that search mixes is then mixed in the column's own state, conserving h_L and q_t per column. 0 = off (2026-11-09, ATSAT row SAT-SPINUP-DRIFT): the column-local test fires single rows on a 1-deg grid; on ATSAT (satchk/mcasmooth_2026-11-09, HA 500 -> 1000) =2 removed the grid-row T steps under the lid (jumps > 3 K at 1000: 164 -> 0) and the 12 K cold rows, MCA still in 100 % of columns. Not measured on this model |
| `ATJUP_THERMAL_WIND` | *follows `ATJUP_NONDIM`* | put the zonal wind's shear into thermal-wind balance at startup |
| `ATJUP_TW_LAT_TAPER` | 10 | latitude band about the equator where the balance is not imposed |
| `ATJUP_TW_MERIDIONAL` | 0 | also adjust v, from the zonal density gradient |
| `ATJUP_HYDRO_SPLIT` | *follows `ATJUP_NONDIM`* | carry the radial buoyancy in p_hydro, leaving p_dyn the rest |
| `ATJUP_BUOY_RHO_LEVEL` | *follows `ATJUP_HYDRO_SPLIT`* | in p_hydro, the density anomaly relative to its level's mean density at the physical factor g L/u_0² (ATNEPT / ATURAN form); 0 = the old r_mix form with the factor of `ATJUP_ND_BUOY` |
| `ATJUP_PBARO_FOLLOW` | 0 | with `ATJUP_HYDRO_SPLIT`: p_baro takes up the change of the buoyancy integral since the start with a time constant of n RK steps; 0 = frozen |
| `ATJUP_PBARO_FOLLOW_ZONAL` | 1 | only the zonal mean of that change, so a local anomaly (a vortex's cold top) keeps its pressure; 0 = cell by cell |
| `ATJUP_VORTEX` | **1** (since 2026-11-20; =0 restores the zonally symmetric start) | the free Great Red Spot: a geostrophically balanced elliptical anticyclone at 22 S / 180 E (`ATJUP_VORTEX_J` 112, `_K` 180) on the balanced background, no solid body (`ATJUP_NO_SEAMOUNT` 1); `ATJUP_V_MASSBAL_ZONAL` follows it. Closed row JUP-FREE-GRS ran it with `ATJUP_HYDRO_SPLIT`, `ATJUP_VORTEX_THERMAL`, `ATJUP_VORTEX_BASE_I=26`, `ATJUP_VORTEX_BOTTOM_I=18`, `ATJUP_PRESS_LINE`, which stay opt-in. Maps: `tools/render.py` in ATSAT, e.g. `--levels 20,28 --zoom=-40,0,150,210` |
| `ATJUP_VORTEX_THERMAL` | 0 | with `ATJUP_VORTEX`, a profile and `ATJUP_HYDRO_SPLIT`: the vortex's decay with height as a cold top (T′) instead of frozen p_baro, so it is advected with the vortex; `ATJUP_VORTEX_BASE_I=26` keeps T′ near −9 K |
| `ATJUP_VORTEX_LOW_I=<il>` | 0 | with `ATJUP_VORTEX_BOTTOM_I=<ib>`: a lens vortex -- nothing (wind, p′, T′) below i = il, the lower ramp from 0 at il to full at ib (0 = the ramp from the floor, as before). Opt-in (row JUP-FREE-GRS, 2026-10-27): the ramp from the floor left its warm anomaly in the slow deep wind as a fixed foot; L14U / L14P (LOW_I 14, BOTTOM_I 22, to 4000) kept the column together in longitude but did not decide the foot / top question, because the vortex barely moved |
| `ATJUP_HYDRO_REF` | 0 | 0 = integrate up from the base (intended), 1 = down from the top |
| `ATJUP_BUOY_SCALE` | 1 | multiplier on top of the buoyancy's physical value |
| `ATJUP_BUOY_RAMP_ITERS` | 0 | ramp the buoyancy in linearly over n iterations |
| `ATJUP_BUOY_PDYN` | 0 | put p_dyn back into the buoyancy density (the old, unstable reading) |
| `ATJUP_PDYN_UNITS` | 1 | read p_dyn as the nondimensional kinematic pressure it is; 0 reads it as bar |
| `ATJUP_PGRAD_SCALE` | 1 | pressure-gradient scale — **1 is correct**, see *Scaling* |
| `ATJUP_POISSON_METRIC` | 1 | the Laplacian's own metric factors; 0 restores the divergence's |
| `ATJUP_PRESS_WALL` | 1 | dp/dn = 0 at the obstacle in the Poisson stencil |
| `ATJUP_PRESS_SWEEPS` | 1 | relaxation sweeps of the pressure equation per call |
| `ATJUP_TURB_CURV` | 1 | spherical curvature terms in the turbulence production |
| `ATJUP_COSTHE_ABS` | 0 | restore the old, non-reversing cos θ (diagnostic only) |
| `ATJUP_NO_SEAMOUNT` | 1 | build no obstacle: a smooth spherical shell, default since 2026-11-01 (no solid cone in the flow); 0 restores the cone (CONE_GEOM 2) |
| `ATJUP_U_INIT_SCALE` | 0 | scale the prescribed **radial** branch of the initial overturning cells; 1.0 (the default until 2026-10-08) peaks at 40 m/s, 0.025 at 1 m/s. The branch is the whole of the model's velocity divergence: at 1.0 the H2O column loses 3.8 % in 224 iterations, at 0 it closes |
| `ATJUP_JET_REF` | **1** (since 2026-11-23, =0 restores the node table; SHARED `planet/JetRef.h`; row GIANTS-ZONAL-WIND) | Rebuilds every latitude row of the prescribed zonal wind w from an observed cloud-level profile instead of the node table (which was read off de Pater & Lissauer Figs. 4.31 / 4.32 and is off on every giant: `literature/zonalwind_compare.png`). Tropopause value per 1-deg planetocentric row from `JetRef.h` (Jupiter: HST/OPAL 2018, Schmider et al. 2024 Table 2; Saturn: Cassini ISS continuum 350-500 mbar, Garcia-Melendo et al. 2011; Uranus: Sromovsky et al. 2015 Table 7; Neptune: Voyager fit of Sromovsky et al. 1993), the model's own vertical shape through `init_v_or_w`, bottom value kept from the node table (0 on ATJUP / ATSAT; ice giants' equator -65 m/s ATURAN, -440 m/s ATNEPT, tapering to 0 at the next node); applied before the u_0 division. Byte checks `satchk/jetref_2026-11-23/check.out` (knob) and `satchk/jetrefdefault_2026-11-23/check.out` (default flip); 1000-it pairs (bottom 0) in `satchk/jetref_2026-11-23`. |
| `ATJUP_RK_SCALAR_SYNC` | 2 | after the saturation adjustment, copy the condensable fields (1) and the temperature (2) into the time-level copies the Runge-Kutta step integrates from, so the adjustment's condensate and latent heat survive the step; 0 (the default until 2026-10-08) discards them every step |
| `ATJUP_SATADJ_NEWTON` | 1 | saturation adjustment: step towards the saturation value with the Newton length 1/(1+G), G the latent-heat feedback, instead of a fixed half step; 0 (the behaviour until 2026-10-09) left half of the adjusted H2O cells and 95 % of the NH3 cells of the first call unconverged after 15 passes — all of them cold, all-ice cells starting far above their target. The knob lives in the shared `SaturationAdjustment.h`; its default is per planet and only ATJUP's is 1 |
| `ATJUP_SATADJ_CONSERVE` | 1 | saturation adjustment: ice above the melting point melts into the liquid and cools the cell by the latent heat of fusion, and an evaporation step takes from the phase that is there; 0 (the behaviour until 2026-10-09) deletes that ice, a sink of water. Shared `SaturationAdjustment.h`, per-planet default, only ATJUP's is 1 |
| `ATJUP_SATADJ_FREEZE` | 0 | since 2026-11-24, shared `SaturationAdjustment.h`: liquid colder than t_00 freezes at entry (liquid -> ice, warmed by the latent heat of fusion, at most up to t_00; needs CONSERVE; gases without a t_00..t_0 band untouched); row GIANTS-PHASE-OUTSIDE: the adjustment only moved the vapour, so liquid carried below t_00 by transport was never frozen (ATSAT, DL 5000 -> 5100: H2O +203 / NH3 +349 g/m2 per 100 it, nothing removing it). Unset byte-identical to before. |
| `ATJUP_CP_H2`, `ATJUP_CP_GAS` | 0 / 0 | since 2026-11-24, shared `CpH2.h`: cp(T) of the set-up adiabat and the dynamics (row GIANTS-ICE-TP). REFUSED on ATJUP (set-up not ported, `cp_lapse_supported()` false): a non-zero value prints a message and stays 0. Defaults on the other three since 2026-11-24 (ATSAT 1 / 0, ATURAN 2 / 1, ATNEPT 1 / 1). |
| `ATJUP_NH3_ICE_ALF` | −0.888 | exponent of the NH₃ ice saturation curve, equal to the liquid's so that the two curves meet at the triple point; −1.2 (the value until 2026-10-09) put the ice saturation pressure 23 to 30 times too low through the NH₃ cloud range |
| `ATJUP_RESTART_STRIDE` | 100 | iterations between the periodic restart dumps `jup_restart_<iter>.bin` (1.18 GB each since 2026-10-31, written only when the state is clean); 0 = none. `checkpoint_save_iter = -1` in the config switches off the one explicit dump only, not these |
| `ATJUP_RESTART_CARRIED` | 1 | since 2026-10-31 (row SAT-RESTART-REPRO): the restart file (magic `"JUP2"`) also holds the 18 derived 3D arrays and `vel_star` that the next iteration reads before it rewrites them (`restart_arrays_carried()`, planet/FileIO_Jup.cpp), and the resume no longer rebuilds (no forced physics block, no `initCh4LocalRho` / `restoreVar`): a resumed run reproduces the straight one byte for byte, whole 3D state, even or odd resume, with `ATJUP_CHROM` too. Before, it differed in every output file. 0 = write the old `"JUP1"` file (814 MB); a JUP1 file still loads and is rebuilt as before, with a warning |
| `ATJUP_BC_RADIUS_POSITIVE` | 1 | keep the species non-negative on the two radial boundary planes where the 2-point extrapolation overshoots below zero (almost all of it `h2o_cloud`): 1 = floor at zero there, the value the end-of-iteration clamp gives anyway; 2 = fall back to the plain copy; 0 = the behaviour until 2026-10-09. Against 0 the state from iteration 1 on is unchanged with 1; the column budget loses two terms of 35 000 g/m² that cancelled, and the iteration-0 files no longer carry negative values |
| `ATJUP_BC_SPECIES_WALL_COPY` | 1 | default since 2026-10-11, 0 restores. Closes the bottom wall for the eleven species: the plane `i = 0` copies the first cell instead of taking the 2-point extrapolation, whose value lets the centred diffusion stencil at `i = 1` draw a third of the interior gradient through the wall. Measured over 1000 iterations: without it the integrated cells gain +0.27 % of H2O while the extrapolated plane falls; with it they change by -0.003 %. `ATJUP_CWB_DIAG` prints plane and integrated cells apart. The top plane and all other fields are untouched |
| `ATJUP_TB_DIAG=<n>` | off | print-only: every n iterations, which stage of the loop and which term of `rhs_t` changed the level-mean temperature (`planet/TbJup.h`, ported from `ATSAT_TB_DIAG` on 2026-10-12) |
| `ATJUP_TB_MASK=<m>` | 0 | print-only, with `ATJUP_TB_DIAG`: 1 = the budget over the obstacle ring only (fluid cells with a solid face neighbour), 2 = over the far field in the obstacle's rows (>= 40 cells in phi from it); each level normalised by its own cells |
| `ATJUP_RES_FLOOR=<s>` | off | print-only (2026-10-17, shared `planet/Reporting.h`): at every residuum printout a second line with the continuity residuum in a geometry whose sin θ is floored at s, the maximum and its cell, then the five parts (du/dr, 2u/r, dv/dθ / r, cot v / r, dw/dφ / (r sin)) at both reported cells in both geometries. The standard line uses the TRUE sine; a model that clamps the polar metric integrates cos θ / s instead of cot θ poleward of asin(s), 23× less at 89° for s = 0.4. ATJUP's own floor is `sinthe_min()` = 0.55. Unset and set leave every output file identical (satchk/rescheck_2026-10-17) |
| `ATJUP_RES_PARTS` | 0 | print-only (shared `planet/Reporting.h`, from ATURAN 2026-10-16): the five parts of the residuum in the column and the rows of the reported cell. Output identical with and without |
| `ATJUP_PDYN_INIT` | 2 | since 2026-10-12: the dynamic pressure starts in balance with the centrifugal term `(v²+w²)/r` of the radial momentum equation, as on the other three models (0 = zero start, as always before). 224 iterations, 0 against 2: radial wind max 0.068 → 0.058 m/s, radial transport of the H2O column +24 → +7 g/m², continuity residuum 0.029 → 0.030, column change unchanged (−0.059 %); the solver does not keep the field (max `p_dyn` 1.4e-3 at the start, 8e-5 at the end) |
| `ATJUP_VTK_STRIDE` | 0 | write the VTK slices every n iterations; 0 = at every `checkpoint`, which also sets the cadence of the printed extrema. `ATJUP_VTK_STRIDE=100` keeps the printout dense and the files few |
| `ATJUP_DT` | 0.005 | nondimensional timestep, 0.005 = 7 s of Jupiter time (default since 2026-10-20; 0.001 = 1.4 s was the default before). Measured stable up to 0.05; the error grows linearly with dt, 0.005 matches 0.001 at equal physical time to 3-4 digits in the extrema |
| `ATJUP_BC_RADIUS_COPY` | 0 | plain copy instead of extrapolation at the radial walls (diagnostic) |
| `ATJUP_BC_MOUNT_COPY` | 0 | plain copy instead of extrapolation for the scalars in the obstacle surface cells (diagnostic; separates the extrapolation overshoot from the advection scheme in the flank-cooling runaway) |
| `ATJUP_T_LIMITER` | 0 | FCT-style clip of the temperature update to the local bounds of the old state — stops `t` leaving the range of its own initial data at the obstacle |
| `ATJUP_CHEM_GATE_ZERO` | 1 | zero the NH3/H2S/NH4SH reaction rates outside the 200–230 K window instead of freezing the last value |
| `ATJUP_CHEM_MOLAR_CONC` | 0 | plain molar concentration `c = rho·w/m` in the NH4SH rate law instead of the `r_mix/sum_c` renormalisation. **Off by default**: dimensionally right, but it collapses the forward rate by ~6 orders and removes the NH4SH cloud entirely — `A`/`T_d` were calibrated against the inflated concentrations and would have to be re-fitted first |
| `ATJUP_NH4SH_CAP` | 0 | NH4SH ceiling; **off, the bound it uses (`r_max`) is a scratch vector left over from the H2S initialiser** |
| `ATJUP_NO_CLAMP` | 0 | disable the zero floor on the species |

**Diagnostics**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_WPROFILE` | 0 | every n iterations, print max **and** area-weighted mean of u,v,w per radial level |
| `ATJUP_PROBE` | — | `"i,j,k"`: term-by-term momentum budget of one cell, all four RK stages |
| `ATJUP_NANCHECK` | 0 | per-iteration census of non-finite cells, by field and index extent |
| `ATJUP_FPE` | 0 | trap floating-point exceptions (use with gdb) |

The first two are what the geometry and obstacle findings below were measured with. `ATJUP_WPROFILE`
pairs the level maximum with the level mean deliberately: a single global extremum cannot tell a
local artifact from a domain-wide source, and that distinction settled the obstacle question.

`ATJUP_RAD_COUPLING` deserves a note. **1.0 is the physically correct value** — the expression is
the exact nondimensional form of dT/dt = Q/(ρ·cₚ) under this model's scaling. But the timestep is
only ~1.4 s of Jupiter time against a radiative relaxation time of ~10⁷ s, so radiative
equilibration needs on the order of 10⁶ iterations. Larger values are a deliberate **acceleration
factor**, not a correction, and they distort the ratio of radiative to advective timescales.
