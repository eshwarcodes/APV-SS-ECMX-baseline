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

This stage hard-codes scratch paths (`E:\...`) that you must replace with
local paths before running.

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

- MATLAB R2023a or later.
- Toolboxes used: Statistics and Machine Learning (sampling/discretize),
  Image Processing (`exportgraphics`), Parallel Computing
  (`parfor` in `Main_load_files*.m` is optional — the loops also run serially).
- Disk: Stage 0 produces ~10–20 GB of intermediate hourly `.mat` files for a
  full year (one APV-SS configuration).

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

## Data Availability Statement (suggested text)

> The MATLAB source code used to produce the baseline APV-SS simulation
> results in this paper is openly available on GitHub at
> `https://github.com/<USER>/<REPO>` under the MIT License. The repository
> contains the radiative pre-compute, coupled crop–thermal–PV solver, and
> single-day lettuce and tomato growth models (toggled via a top-level flag).
> Meteorological forcing (TMY3 Yuma, AZ) and hourly spectral irradiance
> matrices are available from the corresponding author upon reasonable
> request; the pre-computed hourly intermediates (~tens of GB per
> configuration) are not redistributed but can be regenerated from the
> bundled radiative pre-compute scripts.

Replace `<USER>/<REPO>` with the final GitHub path after the repository is
pushed.

## License

MIT — see `LICENSE`.

## Citation

If you use this code, please cite the accompanying paper (Ravishankar et al.,
*Energy Conversion and Management: X*, 2026, in press) and this repository.
