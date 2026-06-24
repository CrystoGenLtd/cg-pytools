# Interactive tools

## `generate_jobscripts.py`

Generates SGE / SLURM job scripts for CrystoGen (CG) + OCC workflows. It
emits the same family of scripts that ship under [`scripts/SGE`](../SGE) and
[`scripts/SLURM`](../SLURM), but with your own cluster paths, core counts and
array sizes filled in instead of the `YOUR_PATH` / `MY_CIF.cif` placeholders.

No dependencies beyond the Python 3 standard library.

> **Note:** the use cases here are currently *hardcoded* — the tool knows a
> fixed set of workflows (the ones under [`scripts/SGE`](../SGE) /
> [`scripts/SLURM`](../SLURM)) and templates them with your paths and resources.
> It does not yet generate arbitrary pipelines. This is expected to become more
> flexible in the future; for now, if your workflow isn't covered, copy the
> closest example script and edit it by hand.

### Pick three things

| Choice       | Options                                                                                                                          |
|--------------|---------------------------------------------------------------------------------------------------------------------------------|
| `--scheduler`| `slurm`, `sge`                                                                                                                   |
| `--mode`     | `parallel` (multicore), `serial`                                                                                                 |
| `--use-case` | `cg_sub`, `solvent_screen_occ`, `solvent_screen_no_occ`, `mixed_solvent_screen`, `growth_rates`, `growth_modifiers_screen`, `size_extract_grow_twice` |

Run `python generate_jobscripts.py --list` for a one-line description of each
use case.

### Usage

List the use cases:

```bash
python generate_jobscripts.py --list
```

Fully on the command line:

```bash
python generate_jobscripts.py \
    --scheduler slurm --mode parallel --use-case solvent_screen_occ \
    --cif paracetamol.cif --array 179 --cores 8 --job-cores 24 \
    --occ /opt/occ/bin/occ --occ-data /opt/occ/share/ \
    --cg-exe /opt/cg/bin/crystogen --cg-key /opt/cg/bin/CrystoGen.key \
    --venv ~/.venv/cg_screen_env --group -o ./my_run
```

Interactive — prompts for anything not supplied (and only for fields the chosen
use case actually needs):

```bash
python generate_jobscripts.py --interactive
```

If you leave out `--scheduler` or `--use-case` and don't pass `--interactive`,
the tool drops into interactive mode automatically.

### What each use case emits

| Use case                  | Files                                              |
|---------------------------|----------------------------------------------------|
| `cg_sub`                  | `cg_folder_sub.sh`                                 |
| `solvent_screen_occ`      | `job.sh`, `screen.sh`, `submit.sh`                |
| `solvent_screen_no_occ`   | `job.sh`, `screen_no_occ.sh`, `submit.sh`         |
| `mixed_solvent_screen`    | `mixed_screen.sh`                                 |
| `growth_rates`            | `growth_rate_screen.sh`                           |
| `growth_modifiers_screen` | `growth_mod_seq.sh`                               |
| `size_extract_grow_twice` | `job.sh`, `screen.sh`, `screen_stg_2.sh`, `submit_all.sh` |

The `submit.sh` / `submit_all.sh` scripts chain the stages with the right
dependency mechanism for the scheduler (`sbatch --dependency=afterok` for SLURM,
`qsub -hold_jid` for SGE).

### Key options

- `--cores` — cores per array task (parallel mode only; serial requests a single core).
- `--job-cores` — cores for the initial OCC lattice-energy job (`job.sh`).
- `--array` — number of array tasks (defaults: 179 for solvent screens, 20 for growth screens).
- `--partition` — SLURM partition (defaults to `multicore` in parallel mode, `serial` in serial mode).
- `--walltime` — wall-clock limit (`#SBATCH -t` / `#$ -l h_rt=`).
- `--group` / `--no-group` — group symmetry-equivalent interactions with `cg_net.py`, or copy net files directly.
- `--venv` / `--no-venv` — python environment to activate before grouping.

The generated scripts are written to `--output-dir` (default: current directory)
and made executable. After generating, supply the run-specific input files
(`input.txt`, `addinput.txt`, `solvents.txt`, structure / net / colour files as
required) alongside them — see the matching `README.txt` under
[`scripts/SGE`](../SGE) / [`scripts/SLURM`](../SLURM) for what each use case
expects.
