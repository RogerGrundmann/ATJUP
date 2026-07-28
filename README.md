# ATJUP

Jupiter atmospheric general circulation model based on the numerical framework of the ATOM climate model.
Solves the 3-D Navier-Stokes equations in a spherical shell extending from the water-cloud level
(p = 10 MPa, T = 350 K) to the upper troposphere (p = 10 kPa, T = 110 K), including full cloud
microphysics and chemistry for the H₂/He atmosphere of Jupiter.

As interested in Computational Fluid Dynamics I looked with my eyes at the statement of planetologist that Jupiters 
Great Red Spot was a storm and thought about a fixed icy body just below the top of the atmosphere which could cause 
a wake flow behind such a body and which could fix it at its position. It looks that the code could provide such a 
different view on Jupiters Great Red Spot.

For all relevant data concerning planet Jupiter and its atmosphere the book Plantetary Sciences 
by Imke de Pater and Jack J. Lissauer was inevitable.

---

## Physics & Numerics

- **Domain:** spherical shell, 41 × 181 × 361 grid points (r × θ × φ), ~2.7 million cells
- **Dynamics:** finite-difference discretisation of the 3-D Navier-Stokes equations in spherical coordinates
- **Time integration:** 4th-order Runge-Kutta (inner loop) with a Poisson pressure solver (outer loop)
- **Parallelism:** OpenMP shared-memory threading
- **Obstacle:** an ellipsoidal solid body (the "SeaMount") installed below the top of the atmosphere —
  the fixed obstruction whose wake is the proposed mechanism for the Great Red Spot, and the model's
  only genuine solid boundary. Its staircase flank is numerically delicate: see *Geometry and
  scaling* before reading features that sit directly on it
- **Thermodynamics:**
  - Temperature initialised as a parabolic pole-to-pole profile (zonally uniform), with the
    tropospheric lapse rate **clamped at a tropopause minimum** (default 110 K), so an isothermal
    stratosphere caps the profile instead of the lapse rate running unbounded to the model top
  - Static pressure follows the polytropic relation below the clamp and the isothermal hydrostatic
    law above it, the two being different physical regimes
  - Boussinesq buoyancy, formed as an anomaly against the area-weighted horizontal mean at each
    level, so only horizontal density contrasts drive vertical motion and the mean is left to
    hydrostatic balance. Its scale is currently inert — see *Geometry and scaling*
  - Clausius-Clapeyron / Sanchez-Lavega SVP formulation for saturation vapour pressures
  - Tao mixed-phase (liquid + ice) saturation adjustment
- **Numerical filtering:** Shapiro de-checkerboarding of the velocity fields, selectable between the
  1-2-1 filter and a 4th-order shear-preserving variant that leaves the resolved zonal-jet shear intact
- **Boundary conditions:** measured temperature/pressure profiles from Voyager (1979), Cassini (2000),
  and HST (2015); rigid radial walls (u = 0) closing the column mass budget for the all-Neumann
  pressure problem; horizontal velocities tapered to a quiet grid ceiling; non-amplifying plain-copy
  treatment at the polar singularity

The three modules below are **off by default**, so a stock run reproduces the baseline dynamics
unchanged. See *Optional modules*.

- **Microphysics:** three-category (rain / snow / graupel) COSMO scheme, run independently for H₂O and
  NH₃, plus Stokes sedimentation of NH₄SH crystals within their stability band. Phase changes at the
  freezing level are formulated in flux space, so melting, freezing and rain evaporation are bounded
  by the hydrometeor flux actually present and cannot create mass or latent heat out of nothing.
  Conversion rates are calibrated against Jupiter's energy budget rather than terrestrial
  microphysical timescales.
- **Radiation:** grey two-stream net-flux solve with H₂/He collision-induced absorption (τ ∝ P²/T),
  additive CH₄/NH₃ bands and a capped cloud/ice continuum. The deep boundary injects Jupiter's
  intrinsic flux (5.4 W/m²); absorbed sunlight enters as a separate shortwave channel — not as a
  downward thermal flux, since the thermal opacity is the wrong absorber for sunlight — with a
  cos(latitude) distribution whose global mean is S(1−A)/4 ≈ 8.3 W/m² against an observed ~8.5.
- **Turbulence:** k-ε (Chien 1982), k-ω (Wilcox 1988) and k-ω SST (Menter 1994), selected by the
  `turb_model` parameter. k and ω are prognostic members of the Runge-Kutta system, so they see the
  same advection, turbulent diffusion and four-stage integration as temperature and the species; the
  eddy viscosity can be fed back into the momentum, heat and species diffusion. The production
  tensor is the full nine-term τ:∇u in the orthonormal spherical basis, curvature terms included.

---

## Chemical Species

| Species | Phases modelled |
|---------|----------------|
| H₂O | vapour · cloud water · cloud ice |
| NH₃ | vapour · cloud · ice |
| CH₄ | vapour · cloud · ice |
| H₂S | vapour |
| NH₄SH | solid crystals (heterogeneous reaction NH₃ + H₂S → NH₄SH at ~230 K) |

---

## Cloud & Circulation Structure

With the microphysics enabled the model reproduces Jupiter's three-deck condensate structure, each
species precipitating where its own saturation curve permits:

| Deck | Species | Form |
|------|---------|------|
| deep, warm | H₂O | rain — liquid, below the freezing level |
| mid-level | NH₄SH | settling crystals, confined to their ~200–230 K stability band |
| aloft | NH₃ | snow only — frozen throughout, its triple point lying below the local temperature |

Zonal jets reach ~130 m/s. With radiation enabled the top-of-atmosphere flux is ~14 W/m² against
Jupiter's observed 13.9.

One point of interpretation: "surface" quantities refer to the **base of the fluid domain**, which
over most of the globe is the deep boundary but over the SeaMount is the obstacle top. Which species
reaches that base therefore depends on where the boundary cuts the vertical structure — over the
obstacle it lies inside the NH₄SH stability band, so the crystals dominate there while H₂O dominates
elsewhere. Read the per-species fields, not only the total.

---

## Repository Layout

```
ATJUP/
├── planet/          # core model (RHS, RK4, thermodynamics, chemistry, radiation,
│                    #   precipitation, turbulence, convective adjustment,
│                    #   boundary conditions, I/O)
├── lib/             # array types, config parser, FFT, utilities
├── cli/             # command-line driver (jup)
├── python/          # Cython bindings (pyatjup)
├── jupiter/         # run directory (XML config, output)
├── tinyxml2/        # vendored XML library
├── param.py         # code-generation script (auto-generates parameter files)
└── Makefile
```

---

## Build

**Dependencies:** C++11 compiler, OpenMP, Python 3, Cython, NumPy.

```bash
make          # generate parameter files, build CLI + Python extension
make jup      # CLI binary only
make python   # Python extension only
make clean
```

The `param.py` script auto-generates several `.inc`/`.pyx`/`.xml` files that parameterise the model;
it runs automatically as part of the build.

Two build behaviours are worth knowing, both of which cause silently stale binaries:

- `python/pyatjup.so` is **not** relinked when only `libatjup.a` changes, so edits under `planet/`
  land in the archive while the running extension stays stale. Force the relink with
  `touch python/pyatjup.cpp && make`.
- `python/pyatjup.cpp` is Cython-generated from `pyatjup.pyx`. After changing `param.py`, **delete**
  it rather than touching it (`rm python/pyatjup.cpp && make python`) — touching makes it newer than
  the `.pyx`, which suppresses re-cythonizing, so new parameters never reach Python. Look for
  `Cythonizing pyatjup.pyx` in the build output to confirm.

---

## Optional modules

Every module is off by default and a stock run is unaffected by their presence. Set the environment
variable to enable.

**Radiation**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_RADIATION` | 0 | grey two-stream solve; fills `Q_rad`, `radiation`, `epsilon` |
| `ATJUP_RAD_COUPLING` | 0 | feed `Q_rad` into the temperature equation (see note below) |
| `ATJUP_SOLAR` | 1 | absorbed shortwave channel (only acts with `ATJUP_RADIATION`) |
| `ATJUP_CIA_STRENGTH` | 1 | scale the H₂/He collision-induced opacity |
| `ATJUP_OPACITY_STRENGTH` | 1 | scale the gas-band and cloud opacity |

**Microphysics**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_PRECIP` | 0 | three-category microphysics + NH₄SH sedimentation |
| `ATJUP_PRECIP_COUPLING` | 0 | feed the latent heat `Q_precip` into the temperature equation |

**Turbulence**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_TURB` | 0 | run the closure |
| `ATJUP_TURB_MODEL` | *param* | override `turb_model` (`k_epsilon`, `k_omega`, `k_omega_SST`) |
| `ATJUP_TURB_COUPLING` | 0 | feed the eddy viscosity into momentum, heat and species diffusion |
| `ATJUP_NUE_MAX` | 1e5 | eddy-viscosity ceiling [m²/s] — a runaway guard, not the operative limiter |
| `ATJUP_ABL_TAPER` | 0 | restore ATOM's boundary-layer taper of the eddy viscosity |

Worth stating plainly, because it is easy to believe otherwise: **`turb_model` in the XML does not
switch the closure on.** It only selects which closure `ATJUP_TURB=1` would run. Without that
variable `tke` and `nue` are identically zero however the XML is written, and without
`ATJUP_TURB_COUPLING=1` the eddy viscosity is computed but never reaches the momentum, heat or
species equations. A run described as "with k-ω SST" that set only the XML had no turbulence in it.

**Initial and boundary conditions**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_IC_TROPO_CLAMP` | 1 | tropopause clamp on the temperature initial condition |
| `ATJUP_T_TROPOPAUSE_MIN` | 110 | tropopause minimum temperature [K] |
| `ATJUP_VEL_SHAPIRO_ORDER` | 2 | 2 = 1-2-1 filter, 4 = shear-preserving |
| `ATJUP_BC_RIGID_LID` | 1 | u = 0 at both radial walls |
| `ATJUP_BC_TOP_TAPER` | 1 | taper v, w to zero over the top three layers |
| `ATJUP_BC_POLE_COPY` | 1 | plain copy instead of extrapolation at the poles |
| `ATJUP_R_NH3_ADD` | *param* | deep well-mixed NH₃ density [kg/m³] |

**Geometry, scaling and the obstacle**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATJUP_METRIC_RADIUS` | 0 (off) | planetary radius [km] for the metric terms — see *Geometry* below |
| `ATJUP_SINTHE_MIN` | 0.55 | polar floor on sin θ; 0.55 means the floor is active poleward of 56.6° |
| `ATJUP_WALL_NUE` | 4 | wall eddy viscosity at the obstacle, in multiples of 1/re; 0 disables |
| `ATJUP_WALL_NUE_LAYERS` | 3 | its ramp depth in cells |
| `ATJUP_LOCAL_RHO` | 0 | use the local density instead of the constant r_mix in q_sat, latent heat and the Lewis groups |
| `ATJUP_NONDIM` | 0 | put the body forces into the units `rhs_u` is written in — see *Scaling* |
| `ATJUP_ND_BUOY` `_COR` `_CENT` | *follow `ATJUP_NONDIM`* | the three factors individually, for attribution |
| `ATJUP_ADIABATIC` | *follows `ATJUP_NONDIM`* | compression work against p_stat, i.e. the dry adiabatic lapse rate |
| `ATJUP_CONV_ADJ` | *follows `ATJUP_NONDIM`* | dry convective adjustment of superadiabatic columns |
| `ATJUP_CONV_ADJ_LAPSE` | 1 | scale its critical lapse rate (1 = the dry adiabat) |
| `ATJUP_CONV_ADJ_PASSES` | 64 | cap on sweeps per column; a warning is printed if it is reached |
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
| `ATJUP_NO_SEAMOUNT` | 0 | build no obstacle at all (diagnostic) |
| `ATJUP_BC_RADIUS_COPY` | 0 | plain copy instead of extrapolation at the radial walls (diagnostic) |
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

---

## Geometry and scaling — read this before interpreting a run

Two structural issues run through the momentum equation. Both are opt-in rather than fixed, because
each is a modelling decision with measured consequences, and neither can be repaired term by term.

**The metric radius.** `rad.z` is built as r₀ = 1.0 with dr = 0.025, so it spans 1.0 to 2.0 across
the shell — while that same dr means one radial step is L_atm/40 = 3.5 km. One number is doing two
incompatible jobs: a vertical grid spacing scaled by L_atm, and the planetary radius that enters
every horizontal metric factor. Jupiter's radius over L_atm is 69911/140 ≈ 499, so every horizontal
derivative is ~500× larger than the geometry warrants and every horizontal Laplacian term 250000×.
`ATJUP_METRIC_RADIUS=69911` gives `rad.z` its geometric meaning; dr is untouched, so the vertical
grid is unchanged and only the curvature moves. The convention is inherited from ATOM, which builds
`rad` identically — this is not specific to ATJUP.

It matters in practice. A fluid cell beside the staircase flank of the SeaMount accelerates without
bound through advective self-amplification, −w ∂w/∂φ growing linearly in w, until the run overflows
around iteration 320. The geometry cuts that driving term by a factor 310 and slows the growth
about sixfold; the wall eddy viscosity (`ATJUP_WALL_NUE`) holds the remainder. Together they give a
run whose maximum vertical velocity *falls* over 500 iterations instead of running away. The
circulation is not otherwise retuned: level means agree within a few percent.

**The body forces were switched off, and `ATJUP_NONDIM=1` switches them on.** `rhs_u` is
nondimensional in units of u₀²/L_atm = 0.0714 m/s², and its body-force terms were never converted
into those units. Each needs its own factor, because each is built from a different combination of
dimensional quantities: buoyancy 1e5·L/u₀² = 1.4e6, Coriolis L/u₀ = 1400, centrifugal L²/u₀² =
1.96e6. `ATJUP_NONDIM=1` applies all three, computed from the configured u₀ and L_atm rather than
written as constants; `ATJUP_ND_BUOY`, `ATJUP_ND_COR` and `ATJUP_ND_CENT` exist for attribution.

**The pressure gradient needs no factor, and the 7.79 that used to be claimed for it was wrong.**
`p_dyn` is not a pressure in bar. Nothing in the model ever assigns it one: it is created solely by
`PressureSolverJup` relaxing ∇²p_dyn = div(aux), and aux is a nondimensional acceleration, so
p_dyn is the nondimensional kinematic pressure p/(ρu₀²) and grad(p_dyn) is already in the units
`rhs_u` wants. That is why raising `ATJUP_PGRAD_SCALE` measured worse rather than better — it was
breaking a balance, not restoring one. Leave it at 1.0. What the factor is genuinely needed for is
the other direction: everywhere p_dyn was **added to p_stat**, which really is in bar, it counted
7.79× too heavily. `p_dyn_to_bar()` now converts it, and the printed field reads 0.002 bar where it
used to read 0.016.

**Three defects came out of doing this properly.**

*The Poisson equation used the metric of a divergence.* Its stencil weights carried 1/r on the θ
term and 1/(r sin θ) on the φ term, where the spherical Laplacian has 1/r² and 1/(r² sin²θ) — the
coefficients of the divergence source twenty lines below, which is where they came from. A
projection's two operators have to be each other's composition. In the original geometry r ≈ 1.5 and
the error is a factor of two, which is how it survived; with `ATJUP_METRIC_RADIUS` the same r is 500.
Corrected (`ATJUP_POISSON_METRIC`, default on), the operator becomes strongly radial, which is the
physical anisotropy for a 140 km shell on a 70000 km planet. Measured over 99 iterations: with the
metric radius it moves p_dyn by 1 % and the velocities not at all; in the original geometry it moves
p_dyn by 30–45 %.

*The buoyancy density contained p_dyn, which is a feedback loop.* With ρ = (p_stat + p_dyn)/(R·T),
buoyancy drives a divergence, the divergence sets p_dyn, p_dyn changes the density and the density
feeds the buoyancy again. In a Boussinesq system the density anomaly is thermodynamic and p_dyn is a
Lagrange multiplier for the velocity constraint; it does not belong in the equation of state. The
loop is invisible while p_dyn is a local smear of the divergence and fatal once the elliptic problem
is actually solved: at full buoyancy, 1 relaxation sweep gives max|u| = 277 m/s after 99 iterations,
50 sweeps 127, and 300 sweeps reaches 69775 m/s by iteration 33 and collapses, with p_dyn at ±60
bar. **More convergence made it worse** — the signature of a feedback, not of a discretisation
error. The buoyancy now reads p_stat alone (`ATJUP_BUOY_PDYN=1` restores the old behaviour), after
which 300 sweeps no longer collapses and simply follows the same convective growth as every other
sweep count — see below.

*The temperature equation had no adiabatic term.* `pressure_t` is the correct nondimensional form of
(1/ρcₚ)Dp/Dt, but it was given p_dyn only. A parcel moving vertically does its work against the
**hydrostatic** pressure, and that term is the entire dry adiabatic lapse rate, dT/dt = −(g/cₚ)w =
2.07 K per kilometre climbed on Jupiter. Without it the model has no static stability at all: a
parcel pushed up keeps its temperature, arrives warmer than its new surroundings and is pushed up
again. This is the reason the buoyancy could not be switched on, and no factor in `rhs_u` would ever
have fixed it — with the buoyancy off, nothing moved vertically for long enough to notice.
`ATJUP_ADIABATIC=1` adds it; at the probe it is 20 % of the temperature tendency.

**What the model does once it can feel all this: it accelerates to about 300 m/s and dies.** A run
with `ATJUP_NONDIM=1 ATJUP_ND_CENT=0 ATJUP_PRESS_SWEEPS=50` grows steadily — max|u| = 80, 129, 171,
241, 323 m/s at iterations 50 to 250 — and goes non-finite at **iteration 233**, in i[10..15], which
is 35 to 52 km. The collapse is not at the pole; the polar maxima that appear at iteration 300 are
already a corrupted field being read.

*Why, and what it is not.* The first explanation tried here was convective: the model does carry a
superadiabatic layer, and it deepens as a run proceeds — two levels near 21–24 km exceeding the dry
adiabat by 0.016 K/km at iteration 100, five levels from 21 to 35 km by up to 0.125 K/km at
iteration 500. That is real and worth knowing. It is also far too weak to be the driver, and a
dry convective adjustment (below) removes it without changing the run at all.

Getting that right depends on one number. **Use the model's own cₚ.** `ChemistryJup` computes
cp_mix = 10655.4 J/(kg·K) for the mixture; the textbook 12000 that looks right for H₂/He moves the
dry adiabat from 2.33 to 2.07 K/km, and 0.26 K/km is twice the whole superadiabatic excess the model
develops. With the wrong value the unstable layer appears to span 14–66 km at 0.39 K/km — three
times too strong and nearly four times too thick.

*What actually drives it is the horizontal contrast.* The buoyancy is an anomaly about the
horizontal mean of each level, so what feeds it is the spread of density **within** a level, and
that spread has never had to be in balance with anything. Measured per level from the reference
run's restart file:

| height | mean T | sd(T) across the level | span(T) | g·sd(ρ)/r_mix | in units of u₀²/L |
|---|---|---|---|---|---|
| 0 km | 318 K | 3.9 K | 21 K | 0.67 m/s² | **9.3** |
| 35 km | 250 K | 13.2 K | 68 K | 0.65 m/s² | **9.1** |
| 77 km | 159 K | 20.8 K | 82 K | 0.63 m/s² | **8.8** |
| 119 km | 111 K | 2.3 K | 18 K | 0.08 m/s² | 1.2 |

A body force of 8–10 against transport terms of order one is not a perturbation. Transport scales as
ũ², so balancing the two needs ũ ≈ 3 in units of u₀ = 100 m/s: **about 300 m/s** — and 323 m/s is
exactly where the run stood when it died. The
model's horizontal temperature structure and its velocity field are simply not solutions of the same
equations, and nothing made them be while the buoyancy was switched off.

Two things follow. The remedy is not another factor: it is to **initialise in thermal-wind balance**,
or to ramp the buoyancy in over the ~2000 iterations geostrophic adjustment needs at this timestep
(`ATJUP_BUOY_RAMP_ITERS`). And the contrast is if anything understated above, because the anomaly is
divided by the constant r_mix = 1.2844 rather than by the level mean density, which at 105 km is
0.083 — proper Boussinesq weighting would make the upper levels 15× stronger still.

The trigger is the temperature going through absolute zero. Tracked every 25 iterations, the
coldest cell in the model holds at −163 °C until iteration 125 and then falls away:

| iteration | 100 | 125 | 150 | 175 | 200 | 225 |
|---|---|---|---|---|---|---|
| min T | −163.2 °C | −163.3 | −163.7 | −164.7 | −168.6 | **−252.9 °C = 20 K** |

Eight iterations later the field is non-finite. A runaway updraft cools adiabatically at 2.33 K per
kilometre climbed, and at several hundred m/s it climbs faster than anything can warm it back. The
alternative candidate — a diffusive limit from a large turbulent viscosity — is ruled out: `tke` and
`nue` are identically zero, because **`ATJUP_TURB` defaults to 0 regardless of the `turb_model` set
in the XML**. With the closure off the only dissipation is 1/re = 0.001 plus the wall viscosity near
the obstacle, so the flow is very nearly inviscid, which is why it accelerates that far.

The practical consequence: **the body forces are correct but the model is not yet ready to run with
them.** `ATJUP_NONDIM` defaults to off, and with it the adiabatic term and the convective
adjustment, so an unchanged command line still gets the previous model.

### The dry convective adjustment

`ConvectiveAdjustmentJup` (Manabe & Strickler 1964) sweeps each column from the bottom, finds every
block steeper than the dry adiabat, and mixes the whole block onto the adiabat while conserving its
mass-weighted enthalpy. It runs after the Runge-Kutta step and its boundary conditions, before the
n-level copies are refreshed, and follows `ATJUP_NONDIM` by default.

Two choices in it are worth knowing. The weights are Δp taken from the hydrostatic `p_stat`, **not a
density** — ρ = p/(R·T) makes ρ·T identically p/R, so a density weight would conserve nothing at all.
And whole segments are mixed rather than adjacent pairs: both converge to the same profile, but
pairwise mixing moves heat one layer per sweep, and a test with `ATJUP_CONV_ADJ_LAPSE=0.5` (which
makes nearly every column qualify) needed 137 million pair operations and still hit the 64-sweep cap
where the segment form settles in two. Conservation is checked in the run itself: the reported
enthalpy drift is 2e-14, which is round-off.

It mixes **temperature only, not composition**, and uses the dry adiabat rather than the saturated
one where cloud is present. Both are deliberate; both are reasonable extensions.

On this model it is close to a no-op, and that is the useful measurement: it fires on 4 to 14 columns
out of 65341 and moves them by 0.03 K, max|u| at iteration 100 is 128.80 with it against 128.78
without, and the run goes non-finite at iteration 245 instead of 233. The superadiabatic layer is
real; it is not what is wrong.

**Timescales, before designing any longer experiment.** The timestep is 1.4 s of Jupiter time, so 500
iterations is twelve minutes. N⁻¹ in the stable layers is 0.08–1 h, and geostrophic adjustment needs
1/f ≈ 2800 s ≈ 2000 iterations. A 99-iteration run covers 139 s and cannot decide a stability
question at all; this is worth remembering against every "stable over N iterations" claim in this
file, including the ones above.

`Forces()` needed two corrections of its own, both in `PresGradForce`: it read p_dyn as bar, and it
divided by `L_atm` where L_atm is in kilometres, making the diagnostic 1000× too large. The other
three components were and remain dimensionally correct — they are force **densities** in N/m³, a
different unit system from the prognostic equation, which is why they carry a 1e5 that `rhs_u` does
not.

**The centrifugal force has a factor but should stay off.** At full strength its radial part is
Ω²R = 2.17 m/s² at the equator, 8.4 % of gravity, and its meridional part peaks near 1 m/s². On the
real planet nothing balances those: they are absorbed into the geopotential, and the answer is
Jupiter's oblateness — the equator sits 4600 km further from the centre than the poles, and the
surfaces of constant effective gravity *are* that shape. This model has a spherical grid, a
spherical lower boundary and a constant radial g, so the meridional part has nothing to work against
and would drive a permanent, entirely spurious pole-to-equator acceleration. The physical treatment
is to fold it into an effective gravity and never write it as a force.

**The polar metric floor.** sin θ is held at `ATJUP_SINTHE_MIN` = 0.55 so that 1/sin and 1/sin²
stay bounded, which means the metric is distorted poleward of 56.6° latitude — 16.5 % of the
sphere. It is a real stability measure, not an oversight: lowering it to 0.15 in the original
geometry produces a non-finite state within 150 iterations, with the entire difference in the polar
caps. With `ATJUP_METRIC_RADIUS` set, the two values become indistinguishable, because the terms the
floor protects are damped by the radius. Lowering the floor is therefore safe only in that
geometry, and 150 iterations is not yet enough evidence to make it the default.

---

## Usage

### Command-line

```bash
./cli/jup jupiter/config_atjup.xml
```

### Python

```python
import sys
sys.path.insert(0, "python")
from pyatjup import Jupiter

model = Jupiter()
model.nm = 100          # iterations
model.checkpoint = 10   # output cadence
model.run()
```

`model.run()` on its own uses the **compiled defaults** and does not read `config_atjup.xml`, so
per-run settings are best made through the exposed properties as above. To read the XML explicitly:

```python
model.load_config(b"config_atjup.xml")   # bytes, not str
```

String-valued parameters likewise read and write as bytes:

```python
model.turb_model = b"k_omega_SST"
```

Output is written as VTK for ParaView: a panorama `.vts` plus radial, meridional and longitudinal
cuts. Alongside the dynamical fields these carry the radiation, precipitation and turbulence
diagnostics, including an all-species surface precipitation map.

Several fields are written in **scaled units** because the VTK writers use fixed-point output with
four decimals, in which the raw SI values would round to zero: the precipitation fluxes `P_*` and
`PrecipSrf_*` in mm/day, `Q_rad` and `Q_precip` in mW/m³, `nue_t` in m²/s, `tke` in m²/s².
The field *names* no longer carry the unit — read this list, not the label.

The condensable species (`h2o`, `nh3`, `ch4`, `h2s`, `nh4sh` and their cloud/ice partners) are mass
**densities** in kg/m³, not mixing ratios. The model had held both readings at once; the density one
is what the initial condition, the microphysics and the sedimentation all assume, and the printout
and the radiative opacity were brought into line with it.

---

## Author

Roger Grundmann — roger.grundmann@web.de
