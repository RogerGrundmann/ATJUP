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
  only genuine solid boundary
- **Thermodynamics:**
  - Temperature initialised as a parabolic pole-to-pole profile (zonally uniform), with the
    tropospheric lapse rate **clamped at a tropopause minimum** (default 110 K), so an isothermal
    stratosphere caps the profile instead of the lapse rate running unbounded to the model top
  - Static pressure follows the polytropic relation below the clamp and the isothermal hydrostatic
    law above it, the two being different physical regimes
  - Boussinesq buoyancy approximation with water-vapour density correction
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
  `turb_model` parameter. k and ω are prognostic, integrated with a positivity-preserving Patankar
  splitting; the eddy viscosity can be fed back into the momentum, heat and species diffusion.
  Turbulent *transport* of k and ω is not yet included, so apart from the SST cross-diffusion term
  they evolve locally.

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

`ATJUP_RAD_COUPLING` deserves a note. **1.0 is the physically correct value** — the expression is
the exact nondimensional form of dT/dt = Q/(ρ·cₚ) under this model's scaling. But the timestep is
only ~1.4 s of Jupiter time against a radiative relaxation time of ~10⁷ s, so radiative
equilibration needs on the order of 10⁶ iterations. Larger values are a deliberate **acceleration
factor**, not a correction, and they distort the ratio of radiative to advective timescales.

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
diagnostics, including an all-species surface precipitation map. Several fields are written in
scaled units with the unit in the field name — `P_rain_mmd` in mm/day, `Q_rad_mW_m3` in mW/m³,
`nue_t_m2s` in m²/s — because the VTK writers use fixed-point output with four decimals, in which
the raw SI values would round to zero.

---

## Author

Roger Grundmann — roger.grundmann@web.de
