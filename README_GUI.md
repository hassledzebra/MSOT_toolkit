# MSOT Pipeline GUI - User Guide

## Overview

The MSOT Pipeline GUI provides an intuitive graphical interface for processing multispectral optoacoustic tomography (MSOT) data. It is designed to replicate the workflow of `examplePipeline.m` while providing an easy-to-use point-and-click interface.

## Features

- **Step-by-step workflow** following the example pipeline
- **Interactive file selection** for input MSOT files and output folders
- **Real-time visualization** of reconstructed and unmixed images
- **Progress monitoring** with status updates and progress bars
- **Flexible configuration** for reconstruction models, resolution, and spectral unmixing
- **Auto-tuning** of speed of sound parameters
- **Multiple endmember selection** for spectral unmixing
- **Results export** to MAT files compatible with downstream analysis

## Requirements

- MATLAB R2021a or later
- Signal Processing Toolbox
- Image Processing Toolbox
- Parallel Computing Toolbox (optional, for parallel processing)

## Quick Start

### Launching the GUI

From MATLAB command window:

```matlab
launchGUI
```

Or directly:

```matlab
MSOTPipelineGUI
```

### Basic Workflow

1. **Setup Tab** - Load your data
   - Click "Browse..." to select your MSOT file (*.msot)
   - Click "Browse..." to select output folder
   - Click "Load Data" to load and view metadata

2. **Reconstruction Tab** - Configure image reconstruction
   - Select Model Type (default: backproject2D)
   - Set image resolution (default: 128x128)
   - Enable auto-tune for speed of sound (recommended)
   - Optionally enable parallel processing

3. **Unmixing Tab** - Select endmembers
   - Check desired endmembers (Hb, HbO2, IR800CW, etc.)
   - Add custom endmembers if needed
   - Select unmixing solver (default: linear)

4. **Process & View Tab** - Run and visualize
   - Click "RUN PROCESSING" to start
   - Watch real-time reconstruction and unmixing
   - Use frame slider to review results
   - Click "Save Results" to export

## Detailed Feature Guide

### Tab 1: Setup

#### Input File Selection
- **MSOT File**: Select the raw MSOT scan file (*.msot format)
- The file contains raw acoustic signals from the scanner

#### Output Folder
- Specify where processed results will be saved
- Folder will be created if it doesn't exist

#### Load Data
- Loads the MSOT file and displays metadata
- Shows scan information including:
  - Number of frames
  - Number of repetitions
  - Sampling frequency
  - Z positions
  - Wavelength information

### Tab 2: Reconstruction Settings

#### Model Type
Choose the reconstruction algorithm:
- **backproject2D** (default): Fast backprojection algorithm, good for quick processing
- **CDMMI**: Coordinate Descent Model-based Method for Iterative reconstruction
- **CDMMI2**: Variant of CDMMI with different parameters
- **dIMMI**: Damped Iterative Model-based Method for Iterative reconstruction

*Note: Iterative methods (CDMMI, dIMMI) provide better image quality but are slower*

#### Image Resolution
- **N_x**: Number of pixels in X direction (default: 128)
- **N_y**: Number of pixels in Y direction (default: 128)
- Range: 64-512 pixels
- Higher resolution = better image detail but slower processing

#### Speed of Sound
The speed of sound in tissue affects image reconstruction accuracy.

- **Auto-Tune** (recommended): Automatically optimizes speed of sound
  - Runs optimization algorithm to find best value
  - Most accurate for your specific data

- **Manual**: Set a specific value (default: 1520 m/s)
  - Use if you know the tissue speed of sound
  - Typical range: 1450-1550 m/s

- **Tune Now**: Manually trigger tuning at any time
  - Uses current resolution settings
  - May take a few moments to compute

#### Parallel Processing
- Enable to use multiple CPU cores
- Requires Parallel Computing Toolbox
- Speeds up processing significantly for large datasets

### Tab 3: Unmixing Settings

Spectral unmixing separates the multispectral data into individual chromophore contributions.

#### Endmember Selection

**Built-in Endmembers:**
- **Deoxyhemoglobin (Hb)**: Deoxygenated blood signal
- **Oxyhemoglobin (HbO2)**: Oxygenated blood signal
- **Water**: Tissue water content
- **Fat**: Lipid signal
- **IR800CW**: Near-infrared fluorescent dye
- **ICG**: Indocyanine green contrast agent

**Custom Endmembers:**
- Enter comma-separated names (e.g., "MyDye1, MyDye2")
- Corresponding spectrum files must exist in `external/spectra/`

#### Unmixing Solver
- **linear**: Standard linear least squares (fastest)
- **nnls**: Non-negative least squares (enforces positive values)
- **nonNegAPCG**: Non-negative accelerated projected conjugate gradient (most accurate)

#### State Filter Type
For temporal smoothing of multispectral data:
- **alphabeta**: Alpha-beta filter (default)
- **alpha**: Alpha filter
- **kalata**: Kalata filter
- **slidingFunction**: Sliding window average

### Tab 4: Process & View

#### Processing Controls

**RUN PROCESSING** button:
- Starts the reconstruction and unmixing pipeline
- Processing cannot be modified once started
- Live progress shown in progress bar and status window

**STOP** button:
- Interrupts processing if needed
- Partial results will be available

#### Progress Monitoring
- **Progress bar**: Visual indicator of completion (0-100%)
- **Frame counter**: Current frame / Total frames
- **Status window**: Timestamped log of processing steps

#### Live Visualization
- **Left panel**: Reconstructed optoacoustic image
  - Shows current frame being processed
  - Color map: jet (blue = low, red = high)

- **Right panel**: Unmixed components
  - Shows separated endmember signals
  - Each component displayed side-by-side

#### Frame Slider
- After processing completes, review any frame
- Drag slider to navigate through frames
- Displays corresponding wavelength information

#### Save Results
- Exports processed data to output folder
- Creates three files:
  1. **msot_processed.mat**: Reconstructed images + wavelengths
  2. **unmix_processed.mat**: Unmixed component images + wavelengths
  3. **pipeline_settings.mat**: Complete processing parameters

## Output Files

### msot_processed.mat
Contains:
- `recon_img`: 4D array (Y × X × Wavelengths × Repetitions)
- `w`: Wavelength array (nm)

### unmix_processed.mat
Contains:
- `unmix_img`: 4D array (Y × X*Components × Wavelengths × Repetitions)
- `w`: Wavelength array (nm)

### pipeline_settings.mat
Contains:
- `pipelineMeta`: Complete structure with all processing parameters
  - `filterSettings`: Signal filtering configuration
  - `reconSettings`: Reconstruction parameters
  - `unmixSettings`: Unmixing configuration
  - `parallelSettings`: Parallel processing settings

## Tips and Best Practices

### For Best Results

1. **Start with defaults**: The default settings work well for most data
2. **Use auto-tune**: Speed of sound auto-tuning improves image quality
3. **Choose appropriate endmembers**: Only select chromophores present in your sample
4. **Check metadata**: Verify scan parameters before processing
5. **Monitor progress**: Watch live images to ensure quality

### Performance Optimization

1. **Lower resolution for testing**: Use 128×128 for quick preview
2. **Enable parallel processing**: Significantly faster for large datasets
3. **Use backproject2D**: Fastest reconstruction method
4. **Selective frame processing**: Process subset of frames first (requires code modification)

### Troubleshooting

**GUI doesn't start:**
- Check MATLAB version (R2021a or later required)
- Verify toolboxes are installed: `ver`
- Check for error messages in command window

**Data loading fails:**
- Verify MSOT file path is correct
- Ensure Java libraries are accessible
- Check file permissions

**Processing errors:**
- Verify at least one endmember is selected
- Check that output folder path is valid
- Ensure sufficient memory for large datasets
- Review status window for specific error messages

**Images look wrong:**
- Try auto-tuning speed of sound
- Verify correct endmembers are selected
- Check that input data quality is good

**Slow processing:**
- Enable parallel processing
- Reduce image resolution
- Use backproject2D model
- Close other MATLAB figures/applications

## Comparison with examplePipeline.m

The GUI implements the same workflow as `examplePipeline.m`:

| examplePipeline.m | GUI Equivalent |
|-------------------|----------------|
| Set file paths manually | Browse buttons in Setup tab |
| `util.MSOTSignalLoader` | "Load Data" button |
| Set `ModelType`, `N_x`, `N_y` | Reconstruction tab settings |
| `recon.backProjectTune` | "Auto-Tune" or "Tune Now" |
| Set `EndmemberNames` | Checkboxes in Unmixing tab |
| Create preFilter, reconstructor, unmixer | Automatic on "Run" |
| For loop with waitbar | Progress bar with live display |
| `save(...)` | "Save Results" button |

## Advanced Usage

### Customizing Processing Parameters

While the GUI exposes the most common parameters, advanced users can modify `MSOTPipelineGUI.m` to add controls for:
- Filter settings (highpass, lowpass frequencies)
- Field of view parameters
- IRF correction settings
- Custom wavelength selection

### Batch Processing

For batch processing multiple scans, consider:
1. Use the GUI to determine optimal parameters
2. Export settings from `pipeline_settings.mat`
3. Write a batch script using `driver.m` or `examplePipeline.m`

### Integration with Analysis Tools

Output MAT files can be loaded in:
- MATLAB for custom analysis
- Python (using scipy.io.loadmat)
- ImageJ/Fiji (export to HDF5 first using existing pipeline tools)

## File Structure

```
MSOT_toolkit/
├── MSOTPipelineGUI.m          # Main GUI application
├── launchGUI.m                # Simple launcher script
├── README_GUI.md              # This documentation
├── examplePipeline.m          # Command-line equivalent
└── +util/, +recon/, +unmix/   # Processing modules
```

## Support and Citation

If you use this GUI in your research, please cite:

> O'Kelly et al., "Microvasculature imaging using model-based reconstruction and spectral unmixing," Scientific Reports, 2021. DOI: 10.1038/s41598-021-97726-1

For issues or questions:
- Review this documentation
- Check MATLAB command window for error messages
- Refer to the original `examplePipeline.m` for workflow comparison

## Version History

### Version 1.0 (2025)
- Initial release
- Tab-based interface
- Real-time visualization
- Auto-tuning support
- Multiple endmember selection
- Progress monitoring
- Results export

## License

This software is provided as-is for research use. See main repository README for license details.
