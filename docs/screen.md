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

The tool ships with `cg-pytools` as the `cg-screen` command, and installing the
package pulls in everything it needs:

```bash
pip install cg-pytools        # or, from a clone: pip install -e .
```

It is implemented in `cgpytools.analysis.screen`, and builds on
`cgpytools.analysis.shape_analysis` (morphology metrics),
`cgpytools.analysis.surfaces` (Wulff / size-file processing),
`cgpytools.io.net` (energy networks), `cgpytools.plot.plot` (styling) and
`cgpytools.io.log` (logging).

---

## Quick Start

### Basic Usage Pattern
```bash
cg-screen -i <path> --<mode> [options]
```

The module is also runnable directly, without the console script:

```bash
python -m cgpytools.analysis.screen -i <path> --<mode> [options]
```

### Minimal Example
```bash
# Analy shapes in general mode
cg-screen -i ./my_shapes --general
```

---

## Analysis Modes

### 1. General Shape Analysis

Analys crystal shapes without solvent-specific context. Best for comparing different crystal morphologies or growth conditions.

#### Command
```bash
cg-screen -i <path> --general [options]
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
cg-screen -i ./crystals --general

# With energy data from CSV
cg-screen -i ./crystals --general --energy_csv energies.csv

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
cg-screen -i <path> --solvent [options]
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
cg-screen -i ./solvent_screen --solvent --solvent_json solvents.json

# With solvent data
cg-screen -i ./solvent_screen --solvent --solvent_json solvents.json --occ --energies

# Exclude specific solvents
cg-screen -i ./solvent_screen --solvent \
    --solvent_json solvents.json --exclude water methanol
```

#### Alternative: Load from Existing Results
```bash
# Re-plot from saved results
cg-screen --results_dir ./RESULTS --solvent
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
cg-screen -i <path> --size --crystallography <json> [options]
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
cg-screen -i ./growth_data --size \
    --crystallography paracetamol_water_cg_results.json

# High-resolution Wulff shapes every 5 steps
cg-screen -i ./growth_data --size \
    --crystallography paracetamol_water_cg_results.json --wulff-interval 10

# With energy data
cg-screen -i ./growth_data --size \
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
cg-screen -i <path> --movies [options]
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
cg-screen -i ./movies --movies

# With labeled trajectories
cg-screen -i ./movies --movies --energy_csv summary.csv
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
cg-screen -i <path> --cda --directions <S> <M> <L> [options]
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
cg-screen -i ./cda_sims --cda --directions " 1  0  0" " 0  1  0" " 0  0  1"

# With energy data
cg-screen -i ./cda_sims --cda \
    --directions " 1  0  0" " 0  1  0" " 0  0  1" --energies
```

#### Alternative: Load from Existing Results
```bash
cg-screen --results_dir ./RESULTS --cda
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
- `--zoom` - Zoom the labelled, `--colour-by` and `--highlight-pic` Zingg plots onto the data.
  By default they keep the full 0-1 Zingg diagram, like the basic Zingg plot, so the shape
  regions stay in proportion
- `--box` - Use bounding box method for Zingg plots (default: SVD)
- `--lmax <value>` - Resolution for Wulff rendering (default: 20, higher = smoother)

### Data Processing
- `--energy_csv <path>` - CSV file with energy data (columns: `shape_name`, `Int_*` or `x*`)
- `--solvent_json <path>` - JSON with solvent properties (for solvent mode)
- `--crystallography <path>` - JSON with crystal structure info (for size mode)
- `--group` - Plot one energy per symmetry group of interactions instead of one per interaction
  line. Each line of a net file is one dimer, and dimers related by symmetry share an ID and so
  an energy, so the per-line `Int_*` columns repeat: in a net file with IDs `1, 2, 2, 3, 3, ...`,
  `Int_2` and `Int_3` are always equal. With `--group` they become `Group_1`, `Group_2`, ...
  numbered by net file ID, on every plot that shows energies, including the correlation
  matrices, `--colour-by energies` and `--highlight`. `group2` or `g2` work as short names for
  `Group_2`. The groups and the lines in each are printed to the log. Applies to `--solvent` and
  `--cda`, whose energies come from net files; the results CSV keeps the `Int_*` columns, so
  `--results_dir` can re-plot either way (the net file is found beside the recorded shape files,
  or under the input folder). Without `--group`, `energies:N` still skips an `Int_*` column that
  repeats one already picked.

  `--group` also exports the **key dimers** if the input folder holds a folder whose name ends
  in `dimers` (e.g. `ce-b3lyp_dimers`) with the dimer XYZ files the net file names
  (`dimer_6.xyz`, ...). They go to `RESULTS/key_dimers/`:
  - `Group_1_dimer_2.xyz`, `Group_2_dimer_1.xyz`, ... - a copy of each group's dimer, named
    after the group and the dimer file it came from
  - `key_dimers.csv` - a lookup table: group, dimer, number of net file lines, centroid
    distance R, the contacts holding the dimer together, the closest intermolecular contact, the
    group's mean/min/max energy across the batch, and where the interaction sits in the crystal
    (below)
  - `key_dimers.png` - every group's dimer drawn side by side at one scale, the two molecules
    in different colours, with contacts marked and described underneath

  If the crystal's CIF is beside the dimers folder (`UREAXX09.cif` next to `UREAXX09_dimers`,
  else a lone CIF there or the first under the input folder), each dimer is also placed in the
  crystal from the symmetry operations on its net file line (`[1A-1A(x,y,z+1)]`):
  - **direction** - the centroid-to-centroid vector as the nearest lattice direction: `along
    [001]` (within 1°), `~[1 1 1] (2.5 deg off)` (within 5°), or else its angles to a/b/c. The
    CSV also gives the fractional vector and the angles to a, b and c.
  - **motif** - what that one interaction builds when every symmetry copy of it is applied: a
    `discrete dimer`, a `closed n-membered ring`, a `chain along [uvw]` (with its repeat
    distance), a `layer parallel to (hkl)` or a `3D network`. A chain is the useful one for
    growth: e.g. urea's double N-H···O interaction is a chain along [001] with a 4.70 Å repeat (c).

  Both go under each picture, which also gets a/b/c axis arrows showing the view (an axis along
  the line of sight is marked ⊙ out of / ⊗ into the page). The arrows are drawn only when the
  dimer file shares the CIF's Cartesian frame, as files written from it do. The molecule each net
  file label (`1A`) stands for is found from the CIF by matching every line's R, so a net file and
  CIF that don't match just leave the directions out, with a warning.

  Contacts are found from the geometry with standard cut-offs, as a quick description rather than
  an energy decomposition: hydrogen bonds (N/O/F-H···N/O/F, H···A ≤ 2.6 Å, angle ≥ 120°,
  counted as single/double/triple...), weak C-H···O/N (≤ 2.7 Å, ≥ 130°), halogen bonds
  (C-Cl/Br/I···A within the vdW sum, ≥ 150°), π stacks (planar 5/6-rings, centroids ≤ 4.4 Å,
  planes within 25°, "offset" if they slip > 1.6 Å) and edge-to-face C-H···π (H ≤ 3.0 Å from a
  ring centroid). A dimer with none of these is described as dispersion only, with its closest
  contact.

### Filtering & Labeling
- `--labels <names>` - Highlight specific items in plots
- `--highlight [columns]` - Automatically label the items with the minimum and maximum of each
  named column on the Zingg plots (`zingg_labeled_<name>.png`, and any `--colour-by` plots),
  e.g. `acetone (max S:M, min dielectric)`. With no columns it uses the aspect ratios (S:M and
  M:L, or the two `AspectRatio_*` columns for CDA). Columns are named as for `--colour-by`, plus
  `ar1`/`ar2` for S:M/M:L, and a `-1` solvent property counts as no value. Each labelled item and
  its value in every named column is printed to the log, so the most extreme shapes can be read
  off without opening the CSV. With two or more columns, each column's min/max pair is drawn in
  its own colour (point, label box and leader line, and the caption on `--highlight-pic`
  pictures), so the pairs can be told apart at a glance; a point that is extreme in several
  columns takes the colour of the first. Combines with `--labels`.
- `--highlight-pic [columns]` - Also save `pictures/zingg_pictures_<name>.png`: the Zingg plot with a
  picture of the crystal behind the minimum and maximum of each named column, e.g.
  `--highlight-pic dielectric acidity`, joined to its point by a leader line. Columns are named as
  for `--highlight` and are independent of it, so the labelled plot and the pictures can show
  different extremes. With no columns it pictures what `--highlight`/`--labels` pick, or the
  aspect ratio extremes if neither is given. Points from `--labels` are always pictured. Extremes
  that differ from `--highlight`'s are printed to the log as "Pictured extremes". Each picture is
  the convex hull of the point's own shape file (the CrystoGen XYZ, or the Wulff mesh for
  `--wulff` and `--size`), drawn from the same angle relative to the crystal's own axes so
  pictures can be compared: longest axis left-right, medium axis up-down, shortest axis towards
  the viewer and tilted to show the thickness. For movies the picture shows the frame in question. At most 12 points are drawn.
  Needs the `shape_file` column the analyses now write to the results CSV, stored relative to
  the input folder; when re-plotting with `--results_dir`, shape files are looked up relative to
  the folder above the results folder (the input folder, with the default `RESULTS` location).
  Not available for CDA, which has no crystal geometry.
- `--colour-by <columns>` (or `--color-by`) - Make one Zingg plot per column, with points
  coloured by that column's value, saved as `coloured_by/zingg_by_<column>_<name>.png`. Any numeric column
  works: solvent properties (`n`, `acidity`, `basicity`, `gamma`, `dielectric`, `aromatic`,
  `halogen`, `solubility`) or interaction energies (`Int_1`, also accepted as `energy1`). Names
  are case-insensitive. Points with no value are shown in grey, and a `-1` solvent property
  (e.g. water's acidity in `solvent.json`) counts as no value. Points from `--labels` and
  `--highlight` are labelled on these plots too. An unknown column is skipped with a warning
  that lists the available ones. Two shortcuts stand for a group of columns: `energies` (one plot
  per interaction energy, `Int_1`, `Int_2`, ...) and `properties` (one per solvent property in the
  results). `energies:N` keeps only N energies, skipping any that are constant. By default
  (`--variance`) these are the N with the largest variance across the batch, most varied first;
  the ones picked and their standard deviations are printed to the log. With `--trend` they are
  instead the N with the strongest trend against the aspect ratios, ranked and logged in the
  same way as `properties:N` below. A large variance shows an energy the solvents change a lot;
  a strong trend shows one that moves with the shape. `properties:N` keeps only the N solvent properties with the
  strongest trend against the aspect ratios: the largest absolute Spearman rank correlation with
  either S:M or M:L, so any steady rise or fall counts, whatever the units. The full ranking, with
  both correlations and the number of solvents behind each, is printed to the log, and
  properties with fewer than 5 values (or one value throughout) are not ranked. A correlation
  from a handful of solvents is a hint, not proof, so check the coloured plot before trusting it.
  The shortcuts work for `--highlight` and `--highlight-pic` too, and
  mix with single columns, e.g. `--colour-by energies:3 dielectric`.
  With `--highlight-pic` as well, the crystal picture plot is also drawn coloured by each
  `--colour-by` column (`pictures/zingg_pictures_by_<column>_<name>.png`), next to the plain one.
- `--overview N` - Quick insight into a system in one flag. Shorthand for
  `--group --highlight energies:N --highlight-pic energies:N --colour-by properties:N energies:N`:
  the N energy groups that vary most across the batch (or trend most with shape, with `--trend`)
  are labelled and their extreme crystals pictured, and those N energies and the N solvent
  properties that trend most with shape each get a coloured Zingg plot and a coloured picture
  plot, plus the plain picture plot and the key dimers. It adds to any of those flags given
  explicitly. Use with `--energies`:
  ```bash
  cg-screen -i . -s solvent.json --solvent --energies --overview 3
  # pick the energies that trend with shape rather than those that vary most
  cg-screen -i . -s solvent.json --solvent --energies --overview 3 --trend
  ```
- `--exclude <names>` - Exclude specific sovlents from analysis

```bash
# Show which solvents give the most extreme shapes from an existing screen
cg-screen --results_dir ./RESULTS --solvent --highlight

# Label the extremes of solvent properties as well as of the aspect ratios
cg-screen --results_dir ./RESULTS --solvent --highlight dielectric acidity ar1 ar2

# Draw the crystals that give the most extreme shapes next to the Zingg plot
cg-screen -i ./solvent_screen --solvent --solvent_json solvents.json --highlight-pic

# ...or the crystals grown in the solvents with the most extreme dielectric and acidity
cg-screen --results_dir ./RESULTS --solvent --highlight-pic dielectric acidity

# ...and how those shapes line up with dielectric constant and the first interaction energy
cg-screen --results_dir ./RESULTS --solvent --highlight --colour-by dielectric energy1

# One Zingg plot per interaction energy and per solvent property
cg-screen --results_dir ./RESULTS --solvent --colour-by energies properties

# ...or just the three interaction energies that differ most between solvents
cg-screen --results_dir ./RESULTS --solvent --colour-by energies:3

# ...or the three solvent properties that best track the shape
cg-screen --results_dir ./RESULTS --solvent --colour-by properties:3

# ...with one energy per symmetry group (Group_1, Group_2, ...) instead of per net file line
cg-screen --results_dir ./RESULTS --solvent --group --colour-by energies:3
```

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
- `coloured_by/` - The `--colour-by` Zingg plots (`zingg_by_<column>_<name>.png`)
- `pictures/` - The `--highlight-pic` crystal picture plots, plain and coloured by each
  `--colour-by` column (`zingg_pictures_<name>.png`, `zingg_pictures_by_<column>_<name>.png`)
- `key_dimers/` - The `--group` key dimer files, table and picture

### Batch Summary
- `batch_summary.log` - What the batch shows, written fresh on every run: the command, then the
  findings printed to the terminal during the run, without the plotting and progress messages.
  That is the energy groups (`--group`), the energies and properties picked by `energies:N` /
  `properties:N` with their standard deviations or the full aspect ratio trend ranking, the
  `--highlight` / `--highlight-pic` extremes with their values, the key dimers with their
  contacts, direction and motif, and the final count of shapes analysed. A short guide at the
  top explains how to read each section. Everything the terminal showed, warnings included, is
  still in `console.log` in the working folder.

---

## Examples

### Example 1: Complete Solvent Screening
```bash
cg-screen \
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
cg-screen \
    -i ./growth_experiment \
    --size \
    --crystallography example_cg_results.json \
    --wulff-interval 10 \
    --energy_csv energies.csv \
```

### Example 3: Multi-Mode Analysis
```bash
# Run general analysis first
cg-screen -i ./study --general --energies

# Then add movie analysis
cg-screen -i ./study --movies

# Finally add CDA
cg-screen -i ./study --cda --directions " 1  0  0" " 0  1  0" " 0  0  1"
```

### Example 4: Re-plotting Existing Results
```bash
# Re-generate plots with different labels
cg-screen \
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
cg-screen -i ./data --general --output_dir ./custom_results
```

### Bounding Box vs SVD Method
```bash
# SVD method (default, more accurate for asymmetric shapes)
cg-screen -i ./data --general

# Bounding box method (faster, good for symmetric shapes)
cg-screen -i ./data --general --box
```

### Consistent Aspect Ratio Scales
```bash
# Force all aspect ratio plots to show 0-1 range
cg-screen -i ./data --general --ar-limits
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
