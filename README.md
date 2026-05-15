# APV-SS Baseline Simulation Code (Yuma, AZ)

MATLAB code accompanying the manuscript on agrivoltaic semi-transparent
(APV-SS) bifacial modules with a wavelength-selective hot mirror, coupled to
lettuce and tomato crop growth and a transient soil-canopy-panel thermal
solver.

This repository contains the **baseline-case** code path for the APV-SS
configuration at the Yuma, AZ climate site. The lettuce and tomato single-day
growth models are both included; a top-level switch selects which crop drives
the coupled simulation.

## Repository layout

```
.
├── run_apvss_baseline.m     Top-level driver. Edit CROP at the top to toggle
│                            between 'lettuce' and 'tomato'.
├── src/                     MATLAB sources (functions + radiative pre-compute
│                            scripts). Added to the MATLAB path by the driver.
├── data/                    Reference inputs (hot-mirror optics, PV module
│                            spectrum, soil + canopy spectra). Large inputs are
│                            documented but NOT bundled — see data/README.md.
├── python/                  Thin Python wrapper (prosail_run.py) that the
│                            MATLAB bridge calls into. Needs upstream pyPro4SAIL
│                            installed separately — see "ProSAIL" below.
├── output/                  Created at runtime. Holds the final baseline
│                            .mat per crop/config.
├── LICENSE                  MIT.
└── README.md                This file.
```

## Pipeline overview

The simulation is run in two stages.

**Stage 0 — radiative pre-compute** (`src/Main_load_files.m`,
`Main_load_files2.m`, `Main_load_files3.m`, `Main_load_files4.m`). For each
sunlit hour of the year, computes hot-mirror angle-resolved optics, ray-traced
ground/panel view factors, spectral irradiance components on the panel front
and rear, and PV operating point (single-diode De Soto model). Outputs are
written to four per-hour `.mat` directories:

| Folder                  | One file per sunlit hour |
|-------------------------|--------------------------|
| `groundIrradiance/`     | spectral PAR / ePAR / NIR on ground grid |
| `panelIrradiance/`      | per-panel front/rear plane-of-array spectra |
| `PVPerformance/`        | per-panel single-diode results, absorbed spectra |
| `viewfactor_groundPV/`  | sparse ground↔panel view factor cells |

Each Stage-0 script opens with an `EDIT THESE PATHS BEFORE RUNNING` config
block (`PLOTS_BASE`, `HOTMIRROR_XLSX`, `SCRATCH_DIR`, and `SCRATCH_DIR1` where
applicable). Defaults are relative (`pwd`-based) — set them to your local
paths or place the files on the MATLAB path before running.
`Main_load_files{2,3,4}.m` are hot-mirror sensitivity variants
(`Hotmirror_T_R_case{2,3,4}.xlsx`); the published baseline uses
`Hotmirror_T_R_original.xlsx` (bundled in `data/`).

**Stage 1 — daily coupling** (`run_apvss_baseline.m`). Loads the hourly
outputs from Stage 0, aggregates to daily totals, builds the active crop
struct (lettuce or tomato), then calls
`simulateCoupledCropThermal4` to step day-by-day through the
Sept 1 → Jan 31 growing window. Each day runs a 3×3 transient thermal solve
(`transientThermalModelSimple2`) coupled to the appropriate single-day crop
model (`simulateLettuceGrowthSingleDay` or
`simulateTomatoGrowthSingleDay_vectorized`), with ET via Penman-Monteith
(`estimateET_mmday`) and PV energy accumulated from the hourly Stage-0 output.

## Crop toggle

`run_apvss_baseline.m`, near the top:

```matlab
CROP = 'lettuce';   % 'lettuce' | 'tomato'
```

Internally a `cropStruct` is built for the selected crop and passed to
`simulateCoupledCropThermal4`, which `switch`-dispatches on `crop.type`. The
lettuce and tomato code paths share the same thermal model and PV rollup;
only the crop-growth state machine differs.

## Requirements

- MATLAB R2023a or later, with a working Python bridge (`pyenv`) — required
  for ProSAIL, see next section.
- Toolboxes used: Statistics and Machine Learning (sampling/discretize),
  Image Processing (`exportgraphics`), Parallel Computing
  (`parfor` in `Main_load_files*.m` is optional — the loops also run serially).
- Disk: Stage 0 produces ~10–20 GB of intermediate hourly `.mat` files for a
  full year (one APV-SS configuration).

## ProSAIL (canopy radiative transfer)

The Stage-0 radiative pre-compute calls the **pyPro4SAIL** package
(H. Nieto; GPL-3) through MATLAB's Python bridge to obtain wavelength-resolved
canopy reflectance and transmittance for ProSAIL leaf parameters fixed in our
`python/prosail_run.py` wrapper.

`python/prosail_run.py` (bundled in this repo) is the **thin wrapper** that
sets the leaf-level Prospect-D inputs used in this paper (N=1.5, Chl=40,
Car=8, EWT=0.01, LMA=0.009, Ant=1.0, hot_spot=0.01, LIDF=(−0.35, −0.15),
default soil) and calls `pypro4sail.four_sail.foursail`. The upstream
**pyPro4SAIL** library itself is **not** redistributed (it is GPL-3 with its
own DOI); install it from source:

```bash
git clone https://github.com/hectornieto/pypro4sail
cd pypro4sail
pip install .
```

Then either:

- set `PROSAIL_PY_DIR` to the directory containing `prosail_run.py`
  (i.e. this repo's `python/` folder, or your own copy), or
- leave it unset — `src/run_prosail_batch.m` defaults to `<repo>/python` when
  `PROSAIL_PY_DIR` is empty.

If you only want to inspect the Stage-1 crop+thermal+PV coupling and skip
the radiative pre-compute, ProSAIL is not needed at runtime; the wrapper
is only invoked from `Main_load_files{,2,3,4}.m`.

Citation:

> Nieto, H. pyPro4SAIL: Vectorized versions of the Prospect5 and 4SAIL
> Radiative Transfer Models. Zenodo, DOI 10.5281/zenodo.11279249.

## Inputs required to run end-to-end

The repository ships function source only. The following inputs are **not
bundled** because of size or because they are platform-specific:

- `plots_base.mat` (~115 MB): meteorology (TMY3), module specs, panel/array
  geometry, crop physiology parameter scalars. Generated locally; see the
  fields loaded in `run_apvss_baseline.m`.
- `Direct_horizn_irradiance_matrix.xlsx` / `Difuse_horizn_irradiance_matrix.xlsx`
  (~60 MB each): hourly spectral DNI/DHI matrices used inside the radiative
  pre-compute. Not required for the Stage-1 driver itself.
- Hourly Stage-0 outputs (`groundIrradiance/`, `panelIrradiance/`,
  `PVPerformance/`, `viewfactor_groundPV/`): produced by Stage-0 scripts.

The small reference files in `data/` (hot-mirror T/R, PV module library,
PVLib spectrum, soil and canopy spectra) **are** bundled so the optical and
electrical functions can be inspected and unit-tested without the full
meteorological dataset.

## How to run

```matlab
>> cd ECMX_paper_code
>> edit run_apvss_baseline.m        % set CROP, CONFIG_NAME, BASE_PATH
>> run_apvss_baseline
```

The driver writes one `baseline_<CONFIG>_<crop>.mat` under `output/`
containing the daily growth state, canopy/soil temperatures,
evapotranspiration, per-panel PV energy, and season totals.

## License

MIT — see `LICENSE`.

## Citation

If you use this code, please cite the accompanying paper (Ravishankar et al.,
*Energy Conversion and Management: X*, 2026, in press) and this repository.
