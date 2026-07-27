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
│                    #   precipitation, turbulence, boundary conditions, I/O)
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
| `ATJUP_BUOY_SCALE` | 1 | buoyancy scale; 1.4e6 is the dimensionally consistent value — see *Scaling* |
| `ATJUP_BUOY_RAMP_ITERS` | 0 | ramp the buoyancy in linearly over n iterations |
| `ATJUP_PGRAD_SCALE` | 1 | pressure-gradient scale; 7.79 is the consistent value |
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

**The body forces are effectively switched off.** `rhs_u` is nondimensional in units of u₀²/L_atm =
0.0714 m/s², but its body-force terms were never converted into those units, each falling short by a
different factor: buoyancy 1.4e6, Coriolis L/u₀ = 1400, centrifugal L/u₀² = 14, pressure gradient
1e5/(r_mix·u₀²) = 7.79. Measured at the probe, the Coriolis contribution is ~1e-4 where consistent
scaling gives ~0.067. The centrifugal force is not small physically — Ω²R = 2.17 m/s² at the
equator, 8.4 % of gravity, the force that gives the planet its oblateness — it is simply not being
felt.

Raising them one at a time does not work, and the measurement is unambiguous: buoyancy at its
consistent value alone drives the radial velocity from 38 to 1535 m/s in 150 iterations, and adding
the pressure gradient's own factor makes it worse, not better. The momentum equation, the Poisson
equation and the projection are one system — the solver builds p_dyn from the unscaled divergence
relation, so scaling only the gradient returns an overshoot rather than a balance. Activating the
buoyancy requires one coherent nondimensionalisation across all three. `Forces()`, by contrast, is
dimensionally correct as it stands: it is a diagnostic force **density** in N/m³, a different unit
system from the prognostic equation, which is why it carries a 1e5 that `rhs_u` does not.

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
