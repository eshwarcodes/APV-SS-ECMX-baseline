# Python bridge for ProSAIL

`prosail_run.py` is the thin wrapper called from
`src/run_prosail_batch.m` via MATLAB's `pyrun`. It pins the leaf-level
Prospect-D inputs used in this paper and returns wavelength-resolved canopy
reflectance, transmittance, and absorptance for a given LAI, sky-light
fraction, and sun geometry.

It depends on the upstream **pyPro4SAIL** package, which is GPL-3 and
**not** redistributed here. Install separately:

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

## Citation

Nieto, H. *pyPro4SAIL: Vectorized versions of the Prospect5 and 4SAIL
Radiative Transfer Models.* Zenodo, DOI
[10.5281/zenodo.11279249](https://doi.org/10.5281/zenodo.11279249).
