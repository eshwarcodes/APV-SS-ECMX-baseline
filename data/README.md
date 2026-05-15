# Reference input data

Small reference files used by the radiative and PV code paths. Bundled here
so the optical / electrical functions can be parsed and exercised without
the full meteorological dataset.

| File                             | Used by                                              |
|----------------------------------|------------------------------------------------------|
| `Hotmirror_T_R_original.xlsx`    | `computeHMOpticsDualFace*`, `Main_load_files.m`      |
| `BifacialPV_ART.xlsx`            | Stack reflectance/transmittance for the PV module    |
| `CEC Modules.csv`                | `read_CEC_CSV_Library` → De Soto single-diode params |
| `PVL_Spectrum.xlsx`              | PVLib reference solar spectrum                       |
| `Schott N-BK7.txt`               | Hot-mirror substrate optical constants               |
| `ProSAIL_WetSoil.txt`            | Soil reflectance for ProSAIL canopy spectra          |
| `lettuce.txt`                    | Lettuce light-response reference                     |
| `relative_irradiance.txt`        | Lettuce-model reference irradiance                   |
| `relative_pfd.txt`               | Lettuce-model reference PFD                          |
| `solarpfd.txt`                   | Solar PFD reference                                  |

## Not bundled (size / availability)

The following must be obtained separately or regenerated locally:

- `plots_base.mat` (~115 MB): aggregated meteorology, geometry, parameters.
  Built from TMY3 weather, the bifacial PV module datasheet, the array CAD,
  and the literature-sourced crop physiology scalars.
- `Direct_horizn_irradiance_matrix.xlsx`, `Difuse_horizn_irradiance_matrix.xlsx`
  (~60 MB each): hourly spectral DNI/DHI matrices used by Stage 0.
- Per-hour `.mat` outputs (~10-20 GB total): produced by
  `src/Main_load_files{,2,3,4}.m`. Re-running Stage 0 regenerates them.

Contact the corresponding author for access to the large inputs.
