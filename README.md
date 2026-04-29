# ATJUP

Jupiter atmospheric general circulation model based on the numerical framework of the ATOM climate model.
Solves the 3-D Navier-Stokes equations in a spherical shell extending from the water-cloud level
(p = 10 MPa, T = 350 K) to the upper troposphere (p = 10 kPa, T = 110 K), including full cloud
microphysics and chemistry for the H₂/He atmosphere of Jupiter.

As interested in Computational Fluid Dynamics I looked with my eyes at the statement of planetologist that Jupiters 
Great Red Spot was a storm and thought about a fixed icy body just below the top of the atmosphere which could cause 
a wake flow behind such a body and which could fix it at its position. It looks that the code could provide such a 
different view on Jupiters Great Red Spot.

---

## Physics & Numerics

- **Domain:** spherical shell, 41 × 181 × 361 grid points (r × θ × φ), ~2.7 million cells
- **Dynamics:** finite-difference discretisation of the 3-D Navier-Stokes equations in spherical coordinates
- **Time integration:** 4th-order Runge-Kutta (inner loop) with a Poisson pressure solver (outer loop)
- **Parallelism:** OpenMP shared-memory threading
- **Thermodynamics:**
  - Temperature initialised as a parabolic pole-to-pole profile (zonally uniform)
  - Boussinesq buoyancy approximation with water-vapour density correction
  - Clausius-Clapeyron / Sanchez-Lavega SVP formulation for saturation vapour pressures
  - Tao mixed-phase (liquid + ice) saturation adjustment
- **Microphysics:** two-category ice scheme adapted from the COSMO weather-forecast model
  (rain/snow precipitation via diagnostic column-equilibrium equations)
- **Radiation:** simplified absorptivity model using water-vapour concentration
- **Boundary conditions:** measured temperature/pressure profiles from Voyager (1979), Cassini (2000), and HST (2015)

---

## Chemical Species

| Species | Phases modelled |
|---------|----------------|
| H₂O | vapour · cloud water · cloud ice |
| NH₃ | vapour · cloud · ice |
| H₂S | vapour |
| NH₄SH | vapour (heterogeneous reaction NH₃ + H₂S → NH₄SH at ~230 K) |

---

## Cloud & Circulation Structure

---

## Repository Layout

```
ATJUP/
├── planet/          # core model (RHS, RK4, thermodynamics, chemistry, I/O)
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
# Generate parameter files and build CLI + Python extension
make

# CLI binary only
make jup

# Python extension only
make python

# Clean
make clean
```

The `param.py` script auto-generates several `.inc`/`.pyx` files that parameterise the model;
it runs automatically as part of the build.

---

## Usage

### Command-line

```bash
./cli/jup jupiter/config_atjup.xml
```

### Python

```python
import sys
sys.path.insert(0, "jupiter")
import pyatjup

model = pyatjup.JupiterModel()
model.load_config("jupiter/config_atjup.xml")
model.run()
```

Output is written as VTK files for visualisation in ParaView.

---

## Author

Roger Grundmann — roger.grundmann@web.de
