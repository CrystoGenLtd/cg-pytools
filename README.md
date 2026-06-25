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

```python
from cgpytools import CrystalShape, ShapeAnalyser

crystal = CrystalShape.from_file("crystal.xyz")

analyser = ShapeAnalyser(zingg_method="svd")
analyser.analyse_crystal(crystal)
metrics = analyser.get_all_frame_metrics()
```

## Modules

| Module | Purpose |
|--------|---------|
| `crystal_io` | Read crystal shapes / frames (`CrystalShape`) |
| `shape_analysis` | Morphology metrics, Zingg classification (`ShapeAnalyser`) |
| `surfaces` | Wulff construction and surface/size-file processing |
| `cg_net` | Interaction-energy network parsing (`CGNet`) |
| `plot` | Shared plotting theme and styling helpers |
| `log` | Logging configuration |

## Scripts

Command-line tools live in [`scripts/python/`](scripts/python):

| Script | Depends on | Description |
|--------|------------|-------------|
| `screen.py` | `cgpytools` | Main crystal-shape analysis tool (general / solvent / size / movie / CDA modes). See [`scripts/python/README.md`](scripts/python/README.md). |
| `growth_kinetics.py` | numpy, pandas, matplotlib, scipy | Standalone: time-evolution plots of size data from simulation subfolders. |
| `growth_rates.py` | numpy, pandas, matplotlib | Standalone: growth-rate-vs-supersaturation summaries from `size.csv` files. |

`screen.py` imports `cgpytools`, so install the package first. The two `growth_*`
scripts are single-file and have no dependency on `cgpytools` — they can be copied
and run on their own.

### HPC job scripts

Example SGE and SLURM job scripts for the CrystoGen + OCC workflows live
under [`scripts/HPC/`](scripts/HPC) organised by use case (solvent screens,
growth rates, etc.). Each folder has a `README.txt` describing the inputs
it needs.

Copy the scripts to a folder containing a CIF or other inputs, depending on
which job type is to be submitted.

## License

MIT — see [LICENSE](LICENSE).
