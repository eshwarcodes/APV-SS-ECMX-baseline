# Python helpers

Two files live here, each used at a different stage of the pipeline:

| File                       | Used by         | Purpose                                          |
|----------------------------|-----------------|--------------------------------------------------|
| `prosail_run.py`           | Stage 0 runtime | ProSAIL canopy reflectance/transmittance via MATLAB `pyrun` |
| `smarts_yuma_driver.ipynb` | Offline (once)  | Pre-compute the hourly spectral irradiance matrices fed to Stage 0 |

Both rely on external upstream packages that are **not** redistributed in
this repo for licensing reasons.

## `prosail_run.py` — ProSAIL bridge

Thin wrapper called from `src/run_prosail_batch.m` via MATLAB's `pyrun`.
Pins the leaf-level Prospect-D inputs used in this paper and returns
wavelength-resolved canopy reflectance, transmittance, and absorptance for
a given LAI, sky-light fraction, and sun geometry.

Depends on the upstream **pyPro4SAIL** package (GPL-3; H. Nieto). Install
separately:

```bash
git clone https://github.com/hectornieto/pypro4sail
cd pypro4sail
pip install .
```

MATLAB's Python bridge must be configured (`pyenv` in MATLAB) to point at
the Python interpreter that has `pypro4sail` installed.

By default, `run_prosail_batch.m` looks for `prosail_run.py` in this folder.
You can override the location by setting the `PROSAIL_PY_DIR` environment
variable to a directory containing `prosail_run.py`.

Citation:

> Nieto, H. *pyPro4SAIL: Vectorized versions of the Prospect5 and 4SAIL
> Radiative Transfer Models.* Zenodo, DOI
> [10.5281/zenodo.11279249](https://doi.org/10.5281/zenodo.11279249).

## `smarts_yuma_driver.ipynb` — SMARTS driver

Reproducible notebook that drives SMARTS 2.9.5 to produce the hourly
wavelength-resolved direct and diffuse irradiance matrices consumed by
Stage 0. Pins the atmospheric configuration used in the paper (Yuma, AZ
site, US Standard Atmosphere, S&F_RURAL aerosol model with TAU5=0.1,
ground albedo type 38, CO2=370 ppm, ASTM-G173 reference spectrum).

Depends on **pySMARTS** (NREL, BSD-3) for the Python interface and on a
local SMARTS install for the engine itself.

```bash
pip install pySMARTS
```

SMARTS itself must be obtained from NREL (license forbids redistribution):

- https://www.nrel.gov/grid/solar-resource/smarts.html

Set the `SMARTSPATH` environment variable (or edit the corresponding cell
in the notebook) to the directory containing `smarts295.exe`.

Citations (required by the SMARTS license):

> Gueymard, C. *Parameterized Transmittance Model for Direct Beam and
> Circumsolar Spectral Irradiance.* Solar Energy 71(5):325–346, 2001.
>
> Gueymard, C. *SMARTS, A Simple Model of the Atmospheric Radiative
> Transfer of Sunshine: Algorithms and Performance Assessment.*
> Professional Paper FSEC-PF-270-95. Florida Solar Energy Center, 1995.
