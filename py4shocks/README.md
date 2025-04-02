# Shock Detection and Tracking Tool

A Python-based tool for the detection and tracking of shock in astrophysical hydrodynamic simulations, particularly focused on PLUTO simulation data.

## Overview

This tool analyzes 3D hydrodynamic simulation data from PLUTO to detect, characterize, and track shock fronts. It includes functionality for:

- Reading and processing VTK/VTR files from PLUTO simulations
- Detecting shock cells using the Rankine-Hugoniot jump conditions
- Calculating Mach numbers for identified shocks
- Tracking the evolution of internal shock and dense cloud material over time.
- Exporting results to VTR files and CSV formats for further analysis

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/username/shock-detection-tracking-tool.git
   cd shock-detection-tracking-tool
   ```

2. Create and activate the Conda environment:
   ```bash
   conda create --name shock_analysis --file requirements.txt
   conda activate shock_analysis
   ```

## Requirements

The following Python packages are required:

- NumPy
- Pandas
- PyVista
- SciPy
- PyEVTK
- YAML

A full list of dependencies is provided in the `requirements.txt` file. You can create a Conda environment using:

```bash
conda create --name <env_name> --file requirements.txt
```

## Configuration

The tool uses a YAML configuration file (`config.yaml`) to set parameters for the analysis. Example key configuration parameters include:

```yaml
route_vtk_folder: /path/to/PLUTO/Test_Problems/folder

calculate_shock_cell:
    activated: True
    output_vtr_shocks_folder: /path/to/VTRs/output/folder
    prefix_name_file: HD_shock

constants:
    gamma: 1.67  # Adiabatic index
    mu: 0.6      # Mean molecular weight

tracking_shocks:
    activated: True
    output_tracking_dir: /path/to/CSV/output/folder
    name_tracking_output: shock_track
    threshold_tracer: 0.95
    threshold_rho: 1.8
    threshold_mach: 3.0
```

## Usage

### Shock Detection

To detect shocks in simulation data:

```bash
python shocks_detection.py
```

This will read the VTK files from the specified input folder, calculate Mach numbers for all detected shocks and also only internal shocks, and save the results as VTR files.

You can also provide command-line arguments to override the configuration file:

```bash
python shocks_detection.py --input_route_vtk_folder /path/to/vtk/files --output_vtr_shocks_folder /path/to/output --prefix_name_file HD_shock --gamma 1.67
```

### Shock Tracking

To track the evolution of the internal shock and dense cloud material over time:

```bash
python shocks_tracking.py
```

This will analyze the detected shocks across multiple time steps and generate a CSV file with tracking data.

You can also provide command-line arguments to override the configuration file:

```bash
python shocks_tracking.py --input_route_vtk_folder /path/to/vtk/files --input_folder_mach /path/to/mach/files --name_mach_file HD_shock --output_tracking_dir /path/to/output --name_tracking_output shock_track --threshold_tracer 0.95 --threshold_rho 3.0 --threshold_mach 2.9
```

## Output Files

The tool generates two main types of output:

1. **VTR Files**: Contains the Mach number data for all and cloud material shocks cells.
2. **CSV Files**: Contains tracking data for the internal shocks and dense cloud material.

## Author

- Wladimir Banda Barragán (wbanda@yachaytech.edu.ec)
- Bryan J. Pinargote (bryanpinargote2000@gmail.com)
- Felix Teutloff

