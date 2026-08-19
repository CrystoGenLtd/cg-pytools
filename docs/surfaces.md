# Wulff Shape & Surface Energy Tool - User Guide

`cg-surfaces` builds Wulff shapes from the surface energies of a crystal: it
reads a CrystoGen/OCC results JSON, converts the per-molecule interaction
energy counts into surface energies, expands the facets by symmetry, and
exports the resulting shape as a mesh, a point cloud or a morphology structure
file.

---

## Table of Contents

- [Overview](#overview)
- [Installation & Requirements](#installation--requirements)
- [Input File](#input-file)
- [Quick Start](#quick-start)
- [Commands](#commands)
  - [info](#1-info)
  - [mesh](#2-mesh)
  - [morphology](#3-morphology)
  - [xyz](#4-xyz)
  - [size](#5-size)
- [Common Options](#common-options)
- [Output Files](#output-files)
- [Python API](#python-api)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

---

## Overview

The tool works from one input — a surface energies JSON — and offers:

- **Surface cut summaries**: every `hkl`, its cut (offset) and its surface
  energy in vacuum and in solvent, plus the lowest-energy cut per unique
  Miller index
- **Wulff shapes** as meshes (`.stl`, `.ply`, `.obj`, …)
- **Morphologies** as CIFs, either the structure annotated with its habit or
  the carved crystallite itself
- **Point clouds** reconstructed from a spherical harmonic expansion of the
  Wulff shape, in the same xyz format the shape-analysis tools read
- **Growth trajectories**: a Wulff shape every *n* time steps of a
  CrystoGen `size.csv`

Every export shares one Wulff construction, so a mesh, a point cloud and a
morphology built from the same energies describe the same crystal.

---

## Installation & Requirements

The tool ships with `cg-pytools` as the `cg-surfaces` command, and installing
the package pulls in everything it needs:

```bash
pip install cg-pytools        # or, from a clone: pip install -e .
```

It is implemented in `cgpytools.analysis.surfaces`, and builds on
[chmpy](https://github.com/peterspackman/chmpy) (Wulff construction, spherical
harmonics, CIF I/O), `trimesh` (mesh export) and `cgpytools.io.log` (logging).

---

## Input File

`cg-surfaces` takes the `*_cg_results.json` file that OCC writes when run with
`--surface-energies <n>`. Without that flag the JSON carries the crystal but no
surface cuts, and only [`info`](#1-info) will work.

### What is read

| JSON block | Used for |
|------------|----------|
| `unit cell` (direct / reciprocal matrices) | Converting Miller indices to Cartesian plane normals |
| `space group` (`symmetry_operations`) | Expanding facets to their symmetry equivalents |
| `asymmetric unit` | Rebuilding the crystal structure, needed by `morphology` |
| `surface_energies.facets` | `hkl`, cut offset, area and interaction energy counts per facet |
| `vacuum` / `solvated` | Per-molecule interaction energies in each context |

Three layouts are accepted, so older result files still load: the facet data is
looked for under `surface_cuts`, then at the top level alongside `vacuum`, and
the crystal block under `surface_cuts.surface_energies.crystal`, `crystal` or
`surface_energies.crystal`.

### Surface energies

Facet surface energies are not stored directly — they are computed from the
interaction energy counts of each cut:

```
E_surface(facet) = 0.16604390679 × ½ × Σ (counts × interaction energies) / area
```

The same facet usually appears several times, once per distinct cut through the
cell (`[1 0 1]` at offsets 0.1782 and 0.3218 below, for instance). Which cut
wins is what [`--reduced`](#common-options) and the second table of `info`
decide.

---

## Quick Start

### Basic Usage Pattern
```bash
cg-surfaces <json> <command> [options]
```

The module is also runnable directly, without the console script:

```bash
python -m cgpytools.analysis.surfaces <json> <command> [options]
```

### Minimal Example
```bash
# What is in this file?
cg-surfaces cg_results.json info

# Wulff shape from the solvated energies, written to ./shapes/wulff.stl
cg-surfaces cg_results.json mesh -o shapes
```

The repository ships a working input at
[`resources/outputs/urea_surface_energies.json`](../resources/outputs/urea_surface_energies.json),
which every example below will run against.

---

## Commands

### 1. `info`

Prints the crystal summary, every surface cut, and the lowest-energy cut per
unique Miller index. The one command that works on a JSON without surface cut
data.

#### Command
```bash
cg-surfaces <json> info [--context {vacuum,solvated}]
```

#### Optional Arguments
- `--context {vacuum,solvated}` - Show only these energies, and use them to pick
  the minimum energy cut. Default: show both contexts, pick on the solvated
  energies.

#### What It Does
1. Reports the number of symmetry operations
2. Reports the number of surface cuts and of unique Miller indices
3. Tables every cut — `hkl`, offset, energy — sorted by energy
4. Tables the lowest-energy cut of each unique `hkl`, with its index into the
   loaded facets

#### Example
```bash
cg-surfaces resources/outputs/urea_surface_energies.json info
```

```
Symmetry operations: 8
Surface cuts: 9 (5 unique Miller indices)
Surfaces from: resources/outputs/urea_surface_energies.json
  hkl         cut  E (vacuum)  E (solvated)
  -----------------------------------------
  [1 1 0]  0.5000      0.1591        0.1010
  [1 0 1]  0.1782      0.3066        0.1282
  [1 0 1]  0.3218      0.3066        0.1282
  [2 0 0]  0.0000      0.2250        0.1428
  ...
Lowest energy cut per unique Miller index (solvated):
  hkl      index     cut  E (solvated)
  ------------------------------------
  [1 1 0]      2  0.5000        0.1010
  [1 0 1]      3  0.1782        0.1282
  ...
```

The `index` column is the position of that cut among the loaded facets, which
is what `reduced_copy` keeps and what a `keep=[...]` selection refers to in the
Python API.

---

### 2. `mesh`

Exports the Wulff shape as a mesh file.

#### Command
```bash
cg-surfaces <json> mesh [-o <dir>] [--name <name>] [--suffix <ext>]
                        [--solvation {vacuum,solvated}] [--reduced]
```

#### Optional Arguments
- `-o, --savedir <dir>` - Output directory (default: `.`)
- `--name <name>` - Base name of the output file (default: `wulff`)
- `--suffix <ext>` - Mesh format suffix, passed to trimesh's exporter
  (default: `.stl`; `.ply`, `.obj`, `.glb` and friends also work)
- `--solvation`, `--reduced` - see [Common Options](#common-options)

#### What It Does
1. Converts the stored interaction energies into surface energies
2. Expands the facets by symmetry, keeping the lowest energy per direction
3. Builds the Wulff construction and triangulates it
4. Writes the mesh and logs its volume

#### Example
```bash
# Every cut in the file
cg-surfaces cg_results.json mesh -o shapes --name wulff

# Only the lowest-energy cut per hkl, from the vacuum energies, as PLY
cg-surfaces cg_results.json mesh -o shapes \
    --name wulff_vacuum --suffix .ply --solvation vacuum --reduced
```

```
Wrote shapes/wulff.stl (volume 0.0088)
```

#### Output Files
- `<name><suffix>` - the Wulff mesh

The mesh is in the units of the surface energies, not Angstrom: the Wulff
construction fixes the shape but not the size. Use `morphology` when a physical
size is needed.

---

### 3. `morphology`

Exports the Wulff shape as a crystal structure file. Two modes, which answer
different questions.

#### Command
```bash
cg-surfaces <json> morphology [--mode {faces,cluster}] [--size <value>]
                              [-o <dir>] [--name <name>] [--suffix <ext>]
                              [--cut-molecules] [--padding <A>]
                              [--solvation {vacuum,solvated}] [--reduced]
```

#### Modes

| Mode | What it writes | Size units | Use when |
|------|----------------|------------|----------|
| `faces` (default) | The structure unchanged — cell, space group, asymmetric unit — plus an `_exptl_crystal_face_*` loop describing the habit | mm (default 0.1) | You want a small, symmetry-preserving file that Mercury and similar viewers draw the morphology from |
| `cluster` | The crystallite carved out of the lattice, in P1 in a box of vacuum | Å (default 50) | The atoms of the particle are the point |

#### Optional Arguments
- `--mode {faces,cluster}` - Which of the two above (default: `faces`)
- `--size <value>` - Largest dimension of the morphology, in mm for `faces` and
  Å for `cluster` (default: 0.1 mm / 50 Å)
- `-o, --savedir <dir>` - Output directory (default: `.`)
- `--name <name>` - Base name of the output file (default: `morphology`)
- `--suffix <ext>` - Structure format (default: `.cif`; `.res` also works)
- `--cut-molecules` - `cluster` only: cut molecules at the surface instead of
  keeping them whole. By default a molecule is kept when its centroid falls
  inside the shape, so the particle is made of whole molecules.
- `--padding <A>` - `cluster` only: vacuum margin around the particle
  (default: 5)
- `--solvation`, `--reduced` - see [Common Options](#common-options)

#### Example
```bash
# Habit recorded on the real structure, 0.1 mm across
cg-surfaces cg_results.json morphology -o shapes --name faces --size 0.1

# The particle itself, 40 A across, molecules kept whole
cg-surfaces cg_results.json morphology -o shapes \
    --mode cluster --name cluster --size 40
```

```
Wrote shapes/faces.cif (12 faces, [0.0856 0.0856 0.1   ] mm)
Wrote shapes/cluster.cif (1568 atoms, [31.06 31.06 33.94] A)
```

In `faces` mode the space group is written both by name
(`_symmetry_space_group_name_H-M`, `_symmetry_Int_Tables_number`,
`_symmetry_cell_setting`) and as the full `_symmetry_equiv_pos_as_xyz` loop; the
H-M symbol keeps chmpy's notation (`P2_1` rather than `P 21`), and no Hall
symbol is written.

Only the faces that survive on the hull are written, so the loop lists the faces
the shape actually shows rather than every symmetry equivalent plane fed in:

```
loop_
_exptl_crystal_face_index_h
_exptl_crystal_face_index_k
_exptl_crystal_face_index_l
_exptl_crystal_face_perp_dist
                   1                    1                    0       0.030251986831
                   1                    0                    1       0.038418849351
                   1                   -1                    0       0.030251986831
                   ...
```

#### Output Files
- `<name><suffix>` - the morphology CIF (or `.res`)

Note that `--size` sets the **largest** dimension; the other two follow from the
shape. A `cluster` smaller than the unit cell fails rather than writing an empty
file, and one larger than about two million atoms is refused — reduce `--size`.

---

### 4. `xyz`

Exports the Wulff shape as a point cloud, reconstructed from its spherical
harmonic expansion. This is the format the shape-analysis side of the package
reads, so it is the bridge between a Wulff shape and `cg-screen` /
`ShapeAnalyser`.

#### Command
```bash
cg-surfaces <json> xyz [-o <dir>] [--name <name>] [--suffix <ext>]
                       [--index <n>] [--l-max <n>]
                       [--solvation {vacuum,solvated}] [--reduced]
```

#### Optional Arguments
- `-o, --savedir <dir>` - Output directory (default: `.`)
- `--name <name>` - Base name of the output file (default: `wulff`)
- `--suffix <ext>` - File extension, **without** a leading dot (default: `txt`)
- `--index <n>` - Index appended to the filename (default: 0), useful when
  writing a trajectory
- `--l-max <n>` - Maximum degree of the spherical harmonic expansion
  (default: 20; higher is smoother and slower)
- `--solvation`, `--reduced` - see [Common Options](#common-options)

#### What It Does
1. Builds the Wulff shape and takes its spherical harmonic transform
2. Reconstructs a point cloud from the coefficients
3. Centres the points on their centroid and scales them onto the unit sphere
4. Writes `{name}_{index}.{suffix}`: point count, comment line, then `x y z`
   triples

#### Example
```bash
cg-surfaces cg_results.json xyz -o shapes --name wulff --index 0 --l-max 20
```

```
Wrote shapes/wulff_0.txt
```

```
1010
shapes/wulff_0.txt
0.09781638416216977  1.3533442754776495e-17  -0.9934055063041708
...
```

#### Output Files
- `<name>_<index>.<suffix>` - the point cloud

Because the points are normalised onto the unit sphere, this export carries the
shape but not the size — aspect ratios and Zingg classes survive, absolute
lengths do not.

---

### 5. `size`

Generates a Wulff shape at intervals through a CrystoGen `size.csv`, turning a
growth simulation into a series of shapes.

#### Command
```bash
cg-surfaces <json> size <size.csv> [<size.csv> ...] [-o <dir>] [-n <steps>]
                        [--format <ext>] [--no-expand-symmetry]
                        [--no-reduce-facets]
```

#### Required Files
- One or more `*size.csv` files. The header names the facets and the first
  column is time:

  ```csv
  time, 1  1  0, 1  0  1, 2  0  0, 1  1  1, 0  0  1
  0,0.10,0.13,0.14,0.18,0.28
  1,0.11,0.13,0.15,0.19,0.29
  ```

#### Optional Arguments
- `-o, --output-dir <dir>` - Directory the combined results CSV is written to
  (default: `.`). Per-file outputs go next to each size file instead.
- `-n, --n-steps <n>` - Generate a shape every *n* time steps (default: 10)
- `--format <ext>` - Mesh format of the generated shapes, no leading dot
  (default: `stl`)
- `--no-expand-symmetry` - Do not expand the facets with the symmetry operations
- `--no-reduce-facets` - Do not reduce the facets to unique crystallographic
  directions

#### What It Does
1. Reads the Miller indices from the header of each size file
2. Every *n* rows, takes that row's energies and builds a Wulff shape
3. Writes each shape to `<size file dir>/shapes/wulff_timestep_%06d.<format>`
4. Summarises the run per file, and across files in the output directory

Existing shape files are skipped rather than overwritten, and a timestep that
fails (a degenerate hull, a row of zeros) is logged and skipped so the rest of
the trajectory still runs.

#### Example
```bash
# A shape per timestep from one file
cg-surfaces cg_results.json size ./run/urea_size.csv -o results -n 1

# Every 10 steps across a whole screen
cg-surfaces cg_results.json size ./screen/*/size.csv -o results -n 10 --format ply
```

```
Found 5 facets: ...
Processing 2 time steps, saving every 1 steps
Saved summary: run/shapes/wulff_generation_summary.csv
Saved detailed data: run/shapes/wulff_detailed_data.json
Saved combined Wulff processing results: results/wulff_processing_results.csv
```

#### Output Files
- `<size dir>/shapes/wulff_timestep_*.<format>` - one mesh per sampled timestep
- `<size dir>/shapes/wulff_generation_summary.csv` - timestep, time value, facet
  and energy counts, output path
- `<size dir>/shapes/wulff_detailed_data.json` - the same, plus the `Int_*`
  energies of each timestep
- `<output dir>/wulff_processing_results.csv` - all files combined, keyed by
  `size_file`

The JSON is still required here — it supplies the lattice and symmetry — but the
energies come entirely from the size files, so `--solvation` and `--reduced` do
not apply to this command.

---

## Common Options

Shared by `mesh`, `morphology` and `xyz`:

- `--solvation {vacuum,solvated}` - Which stored energies to build the shape
  from (default: `solvated`)
- `--reduced` - Keep only the lowest-energy cut of each unique Miller index
  before building the shape, chosen on the `--solvation` energies. Without it,
  every cut in the file contributes and the shape is the inner envelope of all
  of them.

`info` takes `--context` instead, which does the same job for its tables.

The positional `<json>` comes **before** the command, and each command's options
come after it:

```bash
cg-surfaces cg_results.json mesh -o shapes --reduced
#           ^ json          ^ command
```

---

## Output Files

| File | Written by | Contents |
|------|------------|----------|
| `<name>.stl` (or `.ply`, `.obj`, …) | `mesh` | Wulff shape mesh, in energy units |
| `<name>.cif` | `morphology --mode faces` | Structure + `_exptl_crystal_face_*` habit, mm |
| `<name>.cif` | `morphology --mode cluster` | Carved crystallite in P1, Å |
| `<name>_<index>.txt` | `xyz` | Normalised point cloud of the shape |
| `shapes/wulff_timestep_*.stl` | `size` | One shape per sampled timestep |
| `shapes/wulff_generation_summary.csv` | `size` | Per-timestep summary |
| `shapes/wulff_detailed_data.json` | `size` | Per-timestep summary + energies |
| `wulff_processing_results.csv` | `size` | Combined results across size files |
| `console.log`, `report.log` | any | Logging output, written to the working directory |

---

## Python API

Everything the CLI does is available on `CrystalWulff`, which is the better
route when the energies come from somewhere other than the JSON:

```python
from cgpytools.analysis.surfaces import CrystalWulff

wulff = CrystalWulff.from_json("cg_results.json")
print(wulff.n_surfaces, wulff.n_molecules, wulff.n_energies)

mesh = wulff.to_mesh(name="wulff", savedir="shapes", solvation="solvated")
coeffs = wulff.to_mesh(save=False, to_return="sph", l_max=20)
wulff.to_xyz(name="wulff", savedir="shapes", index=0)
faces = wulff.to_morphology(mode="faces", size=0.1, savedir="shapes")
cluster = wulff.to_morphology(mode="cluster", size=50.0, savedir="shapes")

# Lowest-energy cut per unique hkl, as a dict or as a new calculator
reduced = wulff.reduced_surface_energies(context="solvated")
calculator = wulff.reduced_copy(context="solvated")

# Facets and energies from elsewhere
mesh = wulff.calculate_wulff_shape(facets, energies, expand_symmetry=True)
```

Useful pieces beyond the CLI:

- `wulff.wulff(...)` - the underlying `WulffConstruction`, shared by every export
- `wulff_polyhedron(construction, size=...)` - the shape as half spaces, with an
  exact `contains()` point test
- `carve_morphology(crystal, polyhedron, ...)` - carve a crystallite from any
  polyhedron
- `keep=[(1, 1, 0), ...]` on `reduce_hkl`, `get_min_energies`,
  `reduced_surface_energies` and `reduced_copy` - restrict the analysis to a
  subset of Miller indices
- `log_surfaces(..., console=True)` - the `info` tables, at INFO level

A runnable walkthrough lives in
[`resources/examples/surfaces.py`](../resources/examples/surfaces.py):

```bash
python resources/examples/surfaces.py
```

Note that importing `cgpytools.analysis.surfaces` configures logging at import
time, which is why the CLI prints without extra setup and why the two log files
appear in the working directory.

---

## Examples

### Example 1: Inspect, then export the reduced shape
```bash
# Which cuts are in the file, and which one wins per hkl?
cg-surfaces urea_cg_results.json info

# Build the shape from just those winning cuts
cg-surfaces urea_cg_results.json mesh -o shapes --name urea_reduced --reduced
```

### Example 2: Vacuum vs solvated habit
```bash
cg-surfaces urea_cg_results.json morphology \
    -o shapes --name urea_vacuum --solvation vacuum --reduced

cg-surfaces urea_cg_results.json morphology \
    -o shapes --name urea_solvated --solvation solvated --reduced
```

### Example 3: Wulff shape into the shape-analysis pipeline
```bash
# Point cloud in the format cg-screen reads
cg-surfaces urea_cg_results.json xyz -o crystals --name urea --l-max 20

# Then characterise it alongside simulated shapes
cg-screen -i ./crystals --general
```

### Example 4: A growth trajectory as meshes
```bash
cg-surfaces urea_cg_results.json size ./growth/*/size.csv \
    -o results -n 10 --format stl
```

### Example 5: A nanoparticle for further simulation
```bash
cg-surfaces urea_cg_results.json morphology \
    --mode cluster --size 80 --padding 10 -o structures --name urea_80A
```

---

## Troubleshooting

### No Surface Cut Data
**Problem**: `info` prints "No surface cut data in <json>", or another command
raises "No facet data available: this instance was created without surface cut
information"
**Solution**: The JSON has the crystal but no surface energies. Re-run OCC with
`--surface-energies <n>` and use the results file it writes.

### No Crystal Structure
**Problem**: "No crystal structure available: this instance was created without
an asymmetric unit"
**Solution**: `morphology` needs the atoms, which come from the `asymmetric
unit` block of the JSON. `mesh` and `xyz` work without it.

### Cluster Is Empty
**Problem**: "No molecules fall inside the Wulff shape: the requested size is
smaller than the unit cell."
**Solution**: Raise `--size`. In `cluster` mode it is in Angstrom, and it has to
be comfortably larger than a unit cell edge.

### Cluster Too Large
**Problem**: "Morphology would hold more than 2000000 atoms"
**Solution**: Reduce `--size`. Atom count grows with the cube of it, so halving
the size is an eightfold cut.

### Ambiguous Energies
**Problem**: "Ambiguity detected: Please specify 'surface=True' or
'surface=False'."
**Solution**: Python API only — the surface and interaction energy arrays happen
to be the same size for this crystal, so the type cannot be inferred. Pass
`surface=True` or `surface=False` explicitly to `to_mesh`, `to_xyz`,
`to_morphology` or `wulff`.

### Energy / Facet Mismatch
**Problem**: "Number of facets and number of energies received does not match!",
or "Couldn't perform Einstein summation"
**Solution**: The energies must have one entry per loaded facet (surface
energies) or match the `(n_molecules, n_energies)` counts array (interaction
energies). Check `wulff.n_surfaces`, `wulff.n_molecules` and `wulff.n_energies`.
If a `size.csv` header lists more or fewer facets than the row has values, the
extra columns are dropped with a warning.

### Degenerate Shape
**Problem**: "Degenerate Wulff shape: zero extent", or a timestep in `size` mode
logged as a Qhull error
**Solution**: The energies for that shape do not bound a polyhedron — commonly a
row of zeros at the start of a growth run, or a single facet family. Skip those
timesteps, or supply energies for more directions.

### Facets Not Found
**Problem**: "Only N facets found! Found: ... Received (M): ..."
**Solution**: A `keep=[...]` selection named Miller indices that are not in the
file. They must match the loaded `hkl` exactly, including sign; use
`wulff.reduce_hkl()` to see what is available.

### Nothing Printed
**Problem**: The Python API runs but prints nothing
**Solution**: The library logs rather than prints. Most detail is at DEBUG;
pass `console=True` to `log_surfaces` / `log_reduced_surfaces`, or configure
`setup_logging(console_level="DEBUG")` from `cgpytools.io.log`.

---

## See Also

- [`docs/screen.md`](screen.md) - `cg-screen`, the batch shape-analysis tool
- [`README.md`](../README.md) - library overview and Python usage
- [`resources/examples/surfaces.py`](../resources/examples/surfaces.py) -
  runnable example
