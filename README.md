# iThera MSOT MATLAB Toolkit

A comprehensive MATLAB toolkit for processing and analyzing multispectral optoacoustic tomography (MSOT) data acquired using iThera Medical MSOT InVision imaging systems equipped with 512-element transducer arrays.

## Features

- **Graphical User Interface (GUI)** - NEW! Easy-to-use point-and-click interface
- **Signal preprocessing and filtering** - Wiener deconvolution, bandpass filtering
- **Image reconstruction** - Multiple algorithms (backproject2D, CDMMI, dIMMI)
- **Spectral unmixing** - Separate chromophore contributions (Hb, HbO2, dyes, etc.)
- **Auto-tuning** - Automatic speed of sound optimization
- **Parallel processing** - Multi-core support for faster processing
- **HDF5 export** - Compatible with ImageJ and other analysis tools

## Quick Start

### Using the GUI (Recommended for New Users)

Launch the graphical interface:

```matlab
launchGUI
```

Follow the intuitive 4-step workflow:
1. **Setup** - Load MSOT file and view metadata
2. **Reconstruction** - Configure model type, resolution, and speed of sound
3. **Unmixing** - Select endmembers (Hb, HbO2, IR800CW, etc.)
4. **Process & View** - Run processing with live visualization and save results

See [README_GUI.md](README_GUI.md) for detailed GUI documentation.

### Using Command-Line Scripts

For automated processing or batch workflows:

```matlab
% See examplePipeline.m for a complete example
examplePipeline
```

## Requirements

- **MATLAB**: R2021a or later (tested)
- **Required Toolboxes**:
  - Signal Processing Toolbox
  - Image Processing Toolbox
  - Parallel Computing Toolbox (optional, for parallel processing)

## Installation

1. Clone or download this repository
2. Add the toolkit to your MATLAB path:
   ```matlab
   addpath(genpath('/path/to/MSOT_toolkit'))
   ```
3. Launch the GUI or run example scripts

## File Structure

```
MSOT_toolkit/
├── MSOTPipelineGUI.m          # Graphical user interface (NEW!)
├── launchGUI.m                # GUI launcher script
├── README_GUI.md              # Detailed GUI documentation
├── examplePipeline.m          # Example command-line workflow
├── driver.m                   # Main driver function
├── drivePreprocessing.m       # Preprocessing step
├── driveRecons.m              # Reconstruction step
├── driveUnmixing.m            # Unmixing step
├── driveHDF5Writing.m         # HDF5 export step
├── +recon/                    # Reconstruction module
├── +unmix/                    # Unmixing module
├── +filter/                   # Signal filtering module
├── +util/                     # Utility functions
└── external/                  # External libraries and spectra
```

## Processing Pipeline

The toolkit implements a 5-step processing pipeline:

1. **Preprocessing** - Load data, configure settings, tune speed of sound
2. **Reconstruction** - Convert acoustic signals to images
3. **Collation** - Merge parallel processing results
4. **Unmixing** - Separate spectral components
5. **Export** - Save to HDF5 format

The GUI combines steps 1-4 into a streamlined workflow.

## Citation

If you use this toolkit in your research, please cite:

**Primary Citation:**

> Han, Z., MacCuaig, W. M., Gurcan, M. N., Claros-Sorto, J., Garwe, T., Henson, C., ... & McNally, L. R. (2023). Dynamic 2-deoxy-D-glucose-enhanced multispectral optoacoustic tomography for assessing metabolism and vascular hemodynamics of breast cancer. *Photoacoustics*, 32, 100531.

**Original Toolkit:**

> O'Kelly, D., Campbell, J., Gerberich, J.L. et al. A scalable open-source MATLAB toolbox for reconstruction and analysis of multispectral optoacoustic tomography data. *Sci Rep* 11, 19872 (2021). https://doi.org/10.1038/s41598-021-97726-1

## Credits

Modifications are made based on the original toolkit by O'Kelly et al. to allow processing of images acquired using MSOT inVision 512-TF (iThera Medical GmbH, Munich, Germany) equipped with a 512-element, toroidal, ring-shaped transducer array.

GUI developed to provide an accessible interface for the processing pipeline.

## Support

- **GUI Documentation**: See [README_GUI.md](README_GUI.md)
- **Example Scripts**: Review `examplePipeline.m` for command-line usage
- **Issues**: Check error messages and verify requirements are met

## License

This software is provided for research use. Please refer to individual file headers for specific licensing information.


