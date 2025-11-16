# Python MSOT Toolkit

This folder contains a lightweight Python reimplementation of the MATLAB-based MSOT toolkit. The Python code mirrors the original processing stages—signal loading, filtering, reconstruction, and spectral unmixing—while remaining dependency-light and easy to extend.

## Structure

- `msot_toolkit/`
  - `loader.py` – HDF5/NPZ loader with a MATLAB-like API and a synthetic data helper.
  - `filtering.py` – Bandpass pre-filter using SciPy.
  - `reconstruction.py` – Minimal backprojection-inspired reconstructor and a utility to tune the speed of sound.
  - `unmixing.py` – Non-negative least squares spectral unmixing.
  - `pipeline.py` – Orchestrates the pipeline and exposes a demo helper.
- `example_pipeline.py` – Command-line example that runs the full pipeline on a provided file.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r Python/requirements.txt
```

## Usage

### Run the pipeline on your own file

```bash
python Python/example_pipeline.py /path/to/scan.npz Python/output --tune-speed
```

The script saves reconstructions, wavelengths, and (optionally) unmixed volumes as `.npy` files in the output directory and writes the exact settings to `pipeline_settings.json`.

### Synthetic demo

You can explore the pipeline without real data using the built-in demo:

```bash
python - <<'PY'
from msot_toolkit import pipeline
pipeline.demo('Python/output_demo', repetition=2)
print('Demo completed. Outputs saved to Python/output_demo')
PY
```

## Configuration file support

`PipelineSettings.from_json` allows driving the pipeline with a JSON file that mirrors MATLAB defaults. A minimal example looks like:

```json
{
  "input_path": "Python/output_demo/synthetic.npz",
  "output_dir": "Python/out",
  "tune_speed": false,
  "recon": {"image_size": [128, 128], "field_of_view": 0.025},
  "filter": {"low_cut": 100000.0, "high_cut": 8000000.0},
  "unmix": {
    "endmember_names": ["Hb", "HbO2"],
    "spectra": {
      "Hb": [0.9, 0.8, 0.7, 0.6, 0.5],
      "HbO2": [0.6, 0.7, 0.8, 0.9, 1.0]
    }
  }
}
```

Save the JSON file and call:

```bash
python - <<'PY'
from msot_toolkit import PipelineSettings, run_pipeline
settings = PipelineSettings.from_json('config.json')
run_pipeline(settings)
PY
```

## Notes

- The reconstructor is intentionally simple to keep runtime fast; it can be swapped with a physics-accurate implementation without changing the pipeline API.
- HDF5/MSOT support loads the first dataset it finds and falls back to synthetic wavelengths if none are stored.
- The unmixing module expects one reconstruction per wavelength; the loader produces a frame list in that shape by default.
