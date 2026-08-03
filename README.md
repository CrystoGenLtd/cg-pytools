# cg-pytools

A small library for analysing crystal morphology from CrystoGen (formerly CrystalGrower) - crystal growth
simulations: shape characterisation (Zingg ratios, aspect ratios, surface
area / volume), Wulff and surface processing, energy-network parsing, and
plotting utilities. It also ships example SGE / SLURM job scripts for running
CrystoGen + OCC workflows (solvent screens, growth rates, growth modifiers, …)
on HPC clusters.

For a full, GUI-driven workflow, see
[CGAspects](https://github.com/CrystoGenLtd/cgaspects) — CrystoGen's official
(PySide6) data-analysis tool for CrystoGen output. `cg-pytools` is
a stripped-down, lightweight alternative: a scriptable library and command-line
utilities for the same kind of analysis, without the GUI. It is currently
limited in scope, and may in future grow to contain the core of CGAspects.

## Installation

```bash
pip install cg-pytools
```

Or, for local development:

```bash
git clone https://github.com/CrystoGenLtd/cg-pytools.git
cd cg-pytools
pip install -e .
```

## Usage

### Shape analysis

`CrystalShape` reads `.xyz` (multi-frame), `.txt`, and mesh files (`.stl`,
`.ply`, `.glb`); `ShapeAnalyser` turns those coordinates into morphology
metrics (bounding-box or SVD/PCA lengths, Zingg aspect ratios, surface area,
volume, surface-area-to-volume ratio, and a shape class).

```python
from cgpytools import CrystalShape, ShapeAnalyser

crystal = CrystalShape.from_file("crystal.xyz")

analyser = ShapeAnalyser(zingg_method="svd")
analyser.analyse_crystal(crystal)          # all frames; pass frame_idx for one
metrics = analyser.get_all_frame_metrics()  # {frame_idx: ShapeMetrics}

first = metrics[0]
print(first.aspect1, first.aspect2, first.shape, first.surface_area_to_volume_ratio)
```

### Wulff shapes from surface energies

`CrystalWulff` loads the crystallographic information (lattice matrices,
symmetry operations, surface cuts and their energies) from a CrystoGen surface
energies JSON file, expands the facets by symmetry, and builds the Wulff shape.

```python
from cgpytools.analysis.surfaces import CrystalWulff

wulff = CrystalWulff.from_json("cg_results.json")
print(wulff.n_surfaces, wulff.n_molecules)

# Wulff shape from the stored energies, written to disk as a mesh
mesh = wulff.to_mesh(name="wulff", savedir="shapes", suffix=".stl", solvation="solvated")
print(mesh.volume)

# Spherical harmonic coefficients instead of a mesh
coeffs = wulff.to_mesh(solvation="vacuum", save=False, to_return="sph", l_max=20)

# Point cloud (xyz) reconstructed from the harmonic expansion
wulff.to_xyz(name="wulff", savedir="shapes", index=0, l_max=20)

# Lowest-energy cut per unique Miller index
reduced = wulff.reduced_surface_energies(context="solvated")
```

Facets and energies from elsewhere (e.g. a CrystoGen `size` file) can be fed in
directly:

```python
import numpy as np
from cgpytools.analysis.surfaces import (
    CrystalWulff,
    generate_wulff_shapes_from_size_file,
    parse_size_file_headers,
)

wulff = CrystalWulff.from_json("cg_results.json")

facets = parse_size_file_headers("size.csv")     # Miller indices from the headers
energies = np.array([...])                        # one energy per facet
mesh = wulff.calculate_wulff_shape(facets, energies, expand_symmetry=True)

# Or sweep a whole simulation: a shape every n time steps, written to size.csv's
# parent directory under shapes/
results = generate_wulff_shapes_from_size_file(wulff, "size.csv", n_steps=10, file_format="stl")
```

### Energy networks

`CGNet` parses a CrystoGen `.net` interaction-energy file into molecules and
interactions, and can group equivalent interactions or swap the energies out.

```python
from cgpytools import CGNet

net = CGNet("crystal.net")
net.parse()

print(net.n_energies, net.n_unique_energies)
print(net.unique_energies)          # {molecule_label: array of energies}

net.group_net("r")                  # group by distance ("r", "mol_type" or "energy")
net.write("grouped.net")
```

### Plotting and logging

```python
from cgpytools.plot.plot import setup_global_style
from cgpytools.io.log import setup_logging

log = setup_logging(basic_level="DEBUG", console_level="INFO")
setup_global_style("publication")   # "modern", "classic", "minimal", "dark", "publication"
```

`setup_logging` attaches a console handler plus `console.log` and `report.log`
file handlers to the root logger. The library logs rather than prints, so
results such as the `info` summary below only reach the console once logging is
configured. Note that `cgpytools.analysis.surfaces` calls `setup_logging` at
import time (which is why its CLI prints without extra setup, and why importing
it creates the two log files); `cgpytools.io.net` does not, so configure
logging yourself to see its messages.

### Command line

Installing the package provides two commands:

| Command | Module | Purpose |
|---------|--------|---------|
| `cg-surfaces` | `cgpytools.analysis.surfaces` | Wulff shapes and surface energies from an OCC results JSON |
| `cg-screen` | `cgpytools.analysis.screen` | Batch crystal-shape analysis and plotting (general / solvent / size / movie / CDA modes) |

`cg-surfaces` takes the `*_cg_results.json` file that OCC writes when run with
`--surface-energies <n>`; without that flag the JSON has no surface cuts, and
only `info` will work.

#### `cg-surfaces`

```bash
# Summary of the crystal, its surface cuts and the lowest-energy cut per hkl,
# with a column for each of the vacuum and solvated energies
cg-surfaces cg_results.json info

# Restrict both tables to one context, and pick the lowest-energy cuts on it
cg-surfaces cg_results.json info --context vacuum

# Export the Wulff shape as a mesh, using only the lowest-energy cuts
cg-surfaces cg_results.json mesh \
    -o shapes --name wulff --suffix .stl --solvation solvated --reduced

# Export the Wulff shape as an xyz point cloud
cg-surfaces cg_results.json xyz \
    -o shapes --name wulff --index 0 --l-max 20

# Generate shapes from one or more size files, every 10 time steps
cg-surfaces cg_results.json size size.csv \
    -o results -n 10 --format stl
```

#### `cg-screen`

Walks an input directory, discovers the relevant CrystoGen outputs and writes
CSV summaries plus Zingg / heatmap / line plots. See
[`docs/screen.md`](docs/screen.md) for the full guide.

```bash
# Morphology of every shape file under ./crystals
cg-screen -i ./crystals --general

# Solvent screen, including OCC solubilities and net-file energies
cg-screen -i ./solvent_screen --solvent -s solvents.json --occ --energies

# Growth dynamics from size.csv files, a Wulff shape every 10 steps
cg-screen -i ./growth_data --size -c example_cg_results.json --wulff-interval 10
```

`net` remains runnable as a module:

```bash
# Group a net file by interaction distance
python -m cgpytools.io.net -i crystal.net -o grouped.net --group -k r
```

## Modules

| Module | Purpose |
|--------|---------|
| `cgpytools.io.crystal` | Read crystal shapes / frames (`CrystalShape`, `Frame`, `Frames`) |
| `cgpytools.analysis.shape_analysis` | Morphology metrics, Zingg classification (`ShapeAnalyser`) |
| `cgpytools.analysis.surfaces` | Wulff construction and surface/size-file processing (`CrystalWulff`); `cg-surfaces` CLI |
| `cgpytools.analysis.screen` | Batch shape-analysis pipeline and plots (`CrystalShapeAnalysisPipeline`); `cg-screen` CLI |
| `cgpytools.io.net` | Interaction-energy network parsing (`CGNet`) |
| `cgpytools.plot.plot` | Shared plotting theme and styling helpers (`setup_global_style`) |
| `cgpytools.io.log` | Logging configuration (`setup_logging`) |

`CrystalShape`, `ShapeAnalyser` and `CGNet` are re-exported from the top-level
`cgpytools` package; everything else is imported from its module.

## Scripts

Standalone scripts live in [`scripts/python/`](scripts/python):

| Script | Depends on | Description |
|--------|------------|-------------|
| `growth_kinetics.py` | numpy, pandas, matplotlib, scipy | Time-evolution plots of size data from simulation subfolders. |
| `growth_rates.py` | numpy, pandas, matplotlib | Growth-rate-vs-supersaturation summaries from `size.csv` files. |

Both are single-file and have no dependency on `cgpytools` — they can be copied
and run on their own. The screening tool that used to live here is now part of
the package, installed as the `cg-screen` command (see
[`docs/screen.md`](docs/screen.md)).

### HPC job scripts

Example SGE and SLURM job scripts for the CrystoGen + OCC workflows live
under [`scripts/HPC/`](scripts/HPC) organised by use case (solvent screens,
growth rates, etc.). Each folder has a `README.txt` describing the inputs
it needs.

Copy the scripts to a folder containing a CIF or other inputs, depending on
which job type is to be submitted.

## Acknowledgements

The Wulff construction and spherical harmonic tools in `surfaces` are built on
[chmpy](https://github.com/peterspackman/chmpy) by
[Peter Spackman](https://github.com/peterspackman). `chmpy` provides
`WulffConstruction` (the dual-space Wulff construction, its mesh and spherical
harmonic transform) and `reconstruct` (point clouds from spherical harmonic
coefficients).

Parts of `cgpytools.analysis.surfaces` were also written by Peter Spackman, or
adapted from `chmpy` — in particular the symmetry expansion and unique-direction
reduction of Miller indices — and are used here with his permission. Those
sections are marked in the source.

## License

Copyright © 2026 Alvin Jenner W and contributors.

`cg-pytools` is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version — the same licence as
[CGAspects](https://github.com/CrystoGenLtd/cgaspects).

It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU Lesser General Public License for more details.

See [COPYING.LESSER](COPYING.LESSER) (LGPLv3) and [COPYING](COPYING) (GPLv3).

`cg-pytools` depends on [chmpy](https://github.com/peterspackman/chmpy), which
is distributed under the GPL-3.0-or-later licence. chmpy is a dependency rather
than vendored code, but note that a distributed work combining the two is
covered by the GPL.
