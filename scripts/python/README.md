# Crystal Shape Analysis Tool - User Guide

A comprehensive tool for analyzing crystal morphology from computational simulations, supporting multiple analysis modes including solvent screening, growth dynamics, and morphology characterization.

---

## Table of Contents

- [Overview](#overview)
- [Installation & Requirements](#installation--requirements)
- [Quick Start](#quick-start)
- [Analysis Modes](#analysis-modes)
  - [General Shape Analysis](#1-general-shape-analysis)
  - [Solvent Screening Analysis](#2-solvent-screening-analysis)
  - [Si File Analysis](#3-size-file-analysis)
  - [Movie Analysis](#4-movie-analysis)
  - [CDA Analysis](#5-cda-analysis)
- [Common Options](#common-options)
- [File Formats](#file-formats)
- [Output Files](#output-files)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

---

## Overview

This tool analys crystal shapes from computational growth simulations, providing:
- **Morphology characterization** via Zingg plots and aspect ratio analysis
- **Solvent screening** to compare crystal shapes across different solvents
- **Growth dynamics** from timestep data or movie files
- **Energy analysis** from interaction energy data
- **Visualization** with interactive and static plots

---

## Installation & Requirements

### Required Python Packages
```bash
pip install numpy pandas matplotlib seaborn plotly mplcursors tqdm natsort
```

### Required Custom Modules
- `shape_analysis` - Core shape analysis functionality
- `cg_net` - Energy network parsing
- `plot` - Plotting utilities
- `surfaces` - Surface and Wulff shape processing
- `log` - Logging configuration

---

## Quick Start

### Basic Usage Pattern
```bash
python screen.py -i <path> --<mode> [options]
```

### Minimal Example
```bash
# Analy shapes in general mode
python screen.py -i ./my_shapes --general
```

---

## Analysis Modes

### 1. General Shape Analysis

Analys crystal shapes without solvent-specific context. Best for comparing different crystal morphologies or growth conditions.

#### Command
```bash
python screen.py -i <path> --general [options]
```

#### Required Files
- XYZ files matching patterns: `*Aspects.XYZ` or `*CGvisualiser.XYZ`

#### Optional Arguments
- `--energy_csv <path>` - CSV file containing energy data for each shape
- `--labels <names>` - Specific shape names to highlight in plots
- `--energies` - Extract interaction energies from net.txt files

#### What It Does
1. Discovers all XYZ files in the input directory
2. Analys crystal morphology (aspect ratios, Zingg classification)
3. Optionally incorporates energy data from CSV or net files
4. Generates Zingg plots and morphology visualizations
5. Saves results to `general_analysis.csv`

#### Example
```bash
# Basic general analysis
python screen.py -i ./crystals --general

# With energy data from CSV
python screen.py -i ./crystals --general --energy_csv energies.csv

```

#### Output Files
- `general_analysis.csv` - Complete morphology data
- `zingg_plot_general.png` - Zingg classification plot
- `aspect_ratio_*.png` - Various aspect ratio visualizations

---

### 2. Solvent Screening Analysis

Compares crystal morphologies grown in different solvents, ideal for solvent selection studies.

#### Command
```bash
python screen.py -i <path> --solvent [options]
```

#### Required Files
- XYZ files with solvent information in path (e.g., `path/to/solvent_water/*.XYZ`)
- `--solvent_json <path>` - JSON file with solvent properties (defaults to solvent.json in current directory)

#### Optional Arguments
- `--occ` - Include s data (inc. solubility) from `*.stdout` files
- `--energies` - Extract interaction energies from net.txt files
- `--exclude <solvents>` - Exclude specific solvents from analysis
- `--labels <solvents>` - Highlight specific solvents in plots

#### What It Does
1. Groups shapes by solvent (extracted from file paths)
2. Loads solvent properties from JSON
3. Optionally loads solubility data from OCC outputs
4. Analys morphology variations across solvents
5. Creates colored Zingg plots by solvent properties
6. Saves results to `cg_analysis.csv`

#### Example
```bash
# Basic solvent screening
python screen.py -i ./solvent_screen --solvent --solvent_json solvents.json

# With solvent data
python screen.py -i ./solvent_screen --solvent --solvent_json solvents.json --occ --energies

# Exclude specific solvents
python screen.py -i ./solvent_screen --solvent \
    --solvent_json solvents.json --exclude water methanol
```

#### Alternative: Load from Existing Results
```bash
# Re-plot from saved results
python screen.py --results_dir ./RESULTS --solvent
```

#### Output Files
- `cg_analysis.csv` - Morphology data with solvent information
- `zingg_plot_cg.png` - Zingg plot colored by solvent
- `colored_zingg_*.png` - Plots colored by solvent properties
- `labeled_zingg_cg.png` - Plot with specific solvents labeled

---

### 3. Size File Analysis

Analys crystal growth dynamics from `size.csv` files, generating Wulff shapes at different growth stages.

#### Command
```bash
python screen.py -i <path> --size --crystallography <json> [options]
```

#### Required Files
- Size files: `*size.csv`
- `--crystallography <path>` - JSON with crystallographic information (OCC output: "*_cg_results.json")

#### Optional Arguments
- `--wulff-interval <n>` - Generate Wulff shape every n steps (default: 10)
- `--lmax <value>` - Resolution for Wulff rendering (default: 20)
- `--energy_csv <path>` - CSV with energy data
- `--labels <names>` - Shapes to highlight

#### What It Does
1. Processes size.csv files containing growth data over time
2. Generates Wulff shapes at regular intervals
3. Analys morphology evolution during growth
4. Creates time-series visualizations
5. Saves intermediate Wulff shapes and analysis

#### Example
```bash
# Basic si analysis
python screen.py -i ./growth_data --size \
    --crystallography paracetamol_water_cg_results.json

# High-resolution Wulff shapes every 5 steps
python screen.py -i ./growth_data --size \
    --crystallography paracetamol_water_cg_results.json --wulff-interval 10

# With energy data
python screen.py -i ./growth_data --size \
    --crystallography paracetamol_water_cg_results.json --energy_csv energies.csv
```

#### Output Files
- `size_wulff_analysis.csv` - Morphology data over time
- `wulff_shapes/` - Directory with generated Wulff meshes
- `growth_evolution_*.png` - Time-series plots

---

### 4. Movie Analysis

Analys growth dynamics from XYZ movie files containing multiple timesteps.

#### Command
```bash
python screen.py -i <path> --movies [options]
```

#### Required Files
- XYZ movie files (multi-frame XYZ files)

#### Optional Arguments
- `--energy_csv <path>` - CSV with energy data
- `--labels <names>` - Shapes to highlight

#### What It Does
1. Parses multi-frame XYZ files
2. Extracts morphology at each timestep
3. Tracks aspect ratio evolution
4. Creates growth trajectory visualizations
5. Saves results to `movie_analysis.csv`

#### Example
```bash
# Basic movie analysis
python screen.py -i ./movies --movies

# With labeled trajectories
python screen.py -i ./movies --movies --energy_csv summary.csv
```

#### Output Files
- `movie_analysis.csv` - Timestep-resolved morphology data
- `growth_trajectories.png` - Evolution plots
- `movie_*.png` - Individual movie visualizations

---

### 5. CDA Analysis

Analys crystal morphology along specific crystallographic directions, useful for directional growth studies.

#### Command
```bash
python screen.py -i <path> --cda --directions <S> <M> <L> [options]
```

#### Required Files
- CDA simulation files: `*simulation_parameters.txt`

#### Required Arguments
- `--directions <S> <M> <L>` - Three crystallographic directions in order: Short, Medium, Long

#### Optional Arguments
- `--energies` - Extract energies from net files
- `--labels <names>` - Shapes to highlight

#### What It Does
1. Loads CDA simulation data
2. Maps dimensions to specified crystallographic directions
3. Analys anisotropic growth patterns
4. Creates direction-specific Zingg plots
5. Saves results to `cda_analysis.csv`

#### Example
```bash
# CDA analysis along [100], [010], [001]
python screen.py -i ./cda_sims --cda --directions " 1  0  0" " 0  1  0" " 0  0  1"

# With energy data
python screen.py -i ./cda_sims --cda \
    --directions " 1  0  0" " 0  1  0" " 0  0  1" --energies
```

#### Alternative: Load from Existing Results
```bash
python screen.py --results_dir ./RESULTS --cda
```

#### Output Files
- `cda_analysis.csv` - Direction-resolved morphology data
- `zingg_plot_cda.png` - CDA-specific Zingg plot
- `directional_growth_*.png` - Direction analysis plots

---

## Common Options

### Input/Output
- `-i <path>` - Directory containing input files (required for new analysis)
- `--results_dir <path>` - Load from existing results and re-plot
- `--output_dir <path>` - Custom output directory (default: `<input_dir>/RESULTS`)

### Visualization
- `--show` - Display plots interactively (default: save only)
- `--ar-limits` - Force aspect ratio axes to 0-1 range
- `--box` - Use bounding box method for Zingg plots (default: SVD)
- `--lmax <value>` - Resolution for Wulff rendering (default: 20, higher = smoother)

### Data Processing
- `--energy_csv <path>` - CSV file with energy data (columns: `shape_name`, `Int_*` or `x*`)
- `--solvent_json <path>` - JSON with solvent properties (for solvent mode)
- `--crystallography <path>` - JSON with crystal structure info (for size mode)

### Filtering & Labeling
- `--labels <names>` - Highlight specific items in plots
- `--exclude <names>` - Exclude specific sovlents from analysis

---

## File Formats

### Input File Discovery

The tool automatically finds these file types:

| Pattern | Description | Used In |
|---------|-------------|---------|
| `*Aspects.XYZ` | Shape aspect files | General, Solvent, Movies |
| `*CGvisualiser.XYZ` | Visualir files | General, Solvent, Movies |
| `*simulation_parameters.txt` | CDA parameters | CDA |
| `*.*.stdout` | OCC output files | Solvent (with `--occ`) |
| `*size.csv` | Growth size data | Size |

### Energy CSV Format
```csv
shape_name,Int_1,Int_2,x1,x2
shape1,150.5,200.3,0.45,0.62
shape2,175.2,210.1,0.52,0.58
```

### Solvent JSON Format
```json
{
"1-iodopropane": 
    ["1.5058", "0.0000", "0.1500", "41.4500", "6.9626", "0.0000", "0.0000"], 
"benzene": 
    ["1.5011", "0.0000", "0.1400", "40.6200", "2.2706", "1.0000", "0.0000"]
}
```

---

## Output Files

All outputs are saved to the results directory (default: `<input_dir>/RESULTS`).

### CSV Data Files
- `general_analysis.csv` - General mode results
- `cg_analysis.csv` - Solvent screening results
- `cda_analysis.csv` - CDA results
- `size_wulff_analysis.csv` - Size analysis results
- `movie_analysis.csv` - Movie analysis results

### Common Columns
- `name` - Shape/solvent identifier
- `S_length`, `M_length`, `L_length` - Principal axis lengths
- `SR`, `IR` - Short/Intermediate aspect ratios
- `zingg_class` - Zingg classification (prolate, oblate, equant, etc.)
- `volume`, `surface_area` - Geometric properties

### Plot Files
- `zingg_plot_*.png` - Main Zingg classification plots
- `colored_zingg_*.png` - Plots colored by properties
- `labeled_zingg_*.png` - Plots with labeled points
- `aspect_ratio_*.png` - Aspect ratio visualizations
- `*_interactive.html` - Interactive Plotly plots

---

## Examples

### Example 1: Complete Solvent Screening
```bash
python screen.py \
    -i ./solvent_study \
    --solvent \
    --solvent_json solvents.json \
    --occ \
    --energies \
    --labels water ethanol \
    --exclude toluene \
```

### Example 2: Growth Dynamics Study
```bash
python screen.py \
    -i ./growth_experiment \
    --size \
    --crystallography example_cg_results.json \
    --wulff-interval 10 \
    --energy_csv energies.csv \
```

### Example 3: Multi-Mode Analysis
```bash
# Run general analysis first
python screen.py -i ./study --general --energies

# Then add movie analysis
python screen.py -i ./study --movies

# Finally add CDA
python screen.py -i ./study --cda --directions " 1  0  0" " 0  1  0" " 0  0  1"
```

### Example 4: Re-plotting Existing Results
```bash
# Re-generate plots with different labels
python screen.py \
    --results_dir ./RESULTS \
    --solvent \
    --labels methanol octanol
```

---

## Troubleshooting

### No Files Found
**Problem**: "No XYZ files found for general analysis"
**Solution**: 
- Check that files match expected patterns (`*Aspects.XYZ`, `*CGvisualiser.XYZ`)
- Verify `-i` points to correct location
- Ensure files aren't matched by exclusion patterns

### Missing Required Arguments
**Problem**: "Either -i or --results_dir must be specified"
**Solution**: Always provide one of these directory arguments

**Problem**: "--directions must be specified for CDA analysis"
**Solution**: Add `--directions S M L` when using `--cda`

**Problem**: "Can't analyse si files without the crystallographic information"
**Solution**: Provide `--crystallography <json>` when using `--size`

### No Analysis Selected
**Problem**: "No analysis selected"
**Solution**: Add at least one mode flag: `--general`, `--solvent`, `--size`, `--movies`, or `--cda`

### Energy Data Issues
**Problem**: Energy data not loading
**Solution**:
- Check CSV format matches expected columns (`shape_name`, `Int_*` or `x*`)
- Verify shape names match between CSV and XYZ files
- For net files, ensure `net.txt` exists in same directory as shape files

### Solvent Name Extraction
**Problem**: Solvents not detected
**Solution**:
- Ensure directory structure includes solvent name: `path/to/solvent_<name>/shape.XYZ`
- Check solvent names in JSON match extracted names

### Plot Display
**Problem**: Plots not showing
**Solution**: Add `--show` flag to display interactively (default is save-only)

### Performance Issues
**Problem**: Analysis too slow
**Solution**:
- Reduce `--lmax` value (trades quality for speed)
- Increase `--wulff-interval` for si analysis
- Process fewer files at once

---

## Advanced Tips

### Custom Output Location
```bash
python screen.py -i ./data --general --output_dir ./custom_results
```

### Bounding Box vs SVD Method
```bash
# SVD method (default, more accurate for asymmetric shapes)
python screen.py -i ./data --general

# Bounding box method (faster, good for symmetric shapes)
python screen.py -i ./data --general --box
```

### Consistent Aspect Ratio Scales
```bash
# Force all aspect ratio plots to show 0-1 range
python screen.py -i ./data --general --ar-limits
```

---

## Version Info

This guide covers the modular Crystal Shape Analysis Tool with support for:
- General morphology characterization
- Solvent screening
- Growth dynamics (si files and movies)
- CDA directional analysis
- Flexible energy data integration
- Interactive and static visualizations

For questions or issues, consult the tool's logging output or check the generated CSV files for detailed analysis results.
