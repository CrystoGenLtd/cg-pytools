#!/usr/bin/env python3
"""Generate SGE / SLURM job scripts for CrystoGen (CG) + OCC workflows.

This tool emits the same family of scheduler scripts that live under
``scripts/SGE`` and ``scripts/SLURM`` in this repository, but with your own
cluster paths, core counts and array sizes substituted in (instead of the
``YOUR_PATH`` / ``MY_CIF.cif`` placeholders the checked-in examples carry).

Pick a *scheduler* (SGE or SLURM), a *mode* (parallel or serial) and a *use
case*, then either pass everything on the command line or let the tool prompt
you interactively for the bits it still needs.

Examples
--------
List the available use cases::

    python generate_jobscripts.py --list

Fully non-interactive (everything on the command line)::

    python generate_jobscripts.py \
        --scheduler slurm --mode parallel --use-case solvent_screen_occ \
        --cif paracetamol.cif --array 179 --cores 8 --job-cores 24 \
        --occ /opt/occ/bin/occ --occ-data /opt/occ/share/ \
        --cg-exe /opt/cg/bin/crystogen --cg-key /opt/cg/bin/CrystoGen.key \
        --venv ~/.venv/cg_screen_env --group -o ./my_run

Interactive (prompts for anything not supplied)::

    python generate_jobscripts.py --interactive

If you omit ``--scheduler`` or ``--use-case`` and don't pass ``--interactive``,
the tool drops into interactive mode automatically.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# Default placeholder values (mirror the checked-in example scripts).
# ---------------------------------------------------------------------------
DEFAULTS = {
    "cif": "MY_CIF.cif",
    "occ": "YOUR_PATH/occ/bin/occ-0.6.12-linux-x86_64-static/bin/occ",
    "occ_data": "YOUR_PATH/occ/share/",
    "cg_exe": "YOUR_PATH/crystogen-1.2.0-linux-x86_64/bin/crystogen",
    "cg_key": "YOUR_PATH/crystogen-1.2.0-linux-x86_64/bin/CrystoGen.key",
    "venv": "~/.venv/cg_screen_env",
    "cores": 8,
    "job_cores": 12,
    "partition": "multicore",  # slurm partition for parallel mode
    "walltime": "",
}


# ---------------------------------------------------------------------------
# Configuration captured from CLI args / interactive prompts.
# ---------------------------------------------------------------------------
@dataclass
class Config:
    scheduler: str  # "slurm" | "sge"
    mode: str  # "parallel" | "serial"
    use_case: str
    output_dir: Path
    cif: str = DEFAULTS["cif"]
    occ: str = DEFAULTS["occ"]
    occ_data: str = DEFAULTS["occ_data"]
    cg_exe: str = DEFAULTS["cg_exe"]
    cg_key: str = DEFAULTS["cg_key"]
    venv: str = DEFAULTS["venv"]
    cores: int = DEFAULTS["cores"]
    job_cores: int = DEFAULTS["job_cores"]
    array: int = 179
    partition: str = DEFAULTS["partition"]
    walltime: str = DEFAULTS["walltime"]
    group: bool = True


# ---------------------------------------------------------------------------
# Scheduler abstraction.
# ---------------------------------------------------------------------------
@dataclass
class Scheduler:
    key: str
    task_id_var: str
    threads_var: str
    submit_cmd: str


SCHEDULERS = {
    "slurm": Scheduler("slurm", "SLURM_ARRAY_TASK_ID", "SLURM_CPUS_PER_TASK", "sbatch"),
    "sge": Scheduler("sge", "SGE_TASK_ID", "NSLOTS", "qsub"),
}


def header(cfg: Config, sched: Scheduler, *, cores=None, array=None,
           ntasks=None, walltime=None) -> str:
    """Build the scheduler directive header for one script."""
    walltime = walltime or cfg.walltime
    lines = ["#!/bin/bash -l"]
    if sched.key == "slurm":
        lines.append(f"#SBATCH -p {cfg.partition}")
        if ntasks:
            lines.append(f"#SBATCH -n {ntasks}")
        if walltime:
            lines.append(f"#SBATCH -t {walltime}")
        if cores:
            lines.append(f"#SBATCH -c {cores}")
        if array:
            lines.append(f"#SBATCH --array=1-{array}")
    else:  # sge
        lines.append("#$ -cwd")
        if cores:
            lines.append(f"#$ -pe smp.pe {cores}")
        if walltime:
            lines.append(f"#$ -l h_rt={walltime}")
        if array:
            lines.append(f"#$ -t 1-{array}")
    return "\n".join(lines)


def screen_cores(cfg: Config):
    """Per-task cores requested for array/screen jobs.

    Parallel mode requests the chosen core count; serial mode requests a single
    core (SLURM: no ``-c``; SGE: no ``-pe`` line).
    """
    return cfg.cores if cfg.mode == "parallel" else None


# ---------------------------------------------------------------------------
# Reusable shell snippets.
# ---------------------------------------------------------------------------
CIF_BLOCK = """\
# Use the CIF given on the command line, otherwise fall back to a default
if [ -n "$1" ]; then
  CIF=$1
else
  CIF=@@CIF@@
fi"""

CG_EXE_BLOCK = """\
# CrystoGen key + executable
CG_KEY="@@CG_KEY@@"
CG_EXE="@@CG_EXE@@\""""


def venv_block(cfg: Config) -> str:
    if not cfg.venv:
        return ""
    return (
        "\n# Activate the python environment (needed for cg_net.py grouping)\n"
        "source @@VENV@@/bin/activate\n"
    )


def net_block(cfg: Config) -> str:
    """The line(s) that produce ``net.txt`` for a solvent folder."""
    if cfg.group:
        return (
            "# Build net.txt for CG, grouping symmetry-equivalent interactions\n"
            'python cg_net.py -i "$simulation_path/${CIF%.cif}_${SOLVENT}_net.txt"'
            ' -o "$solvent_folder/net.txt" --group'
        )
    return (
        "# Copy the per-solvent net file produced by OCC (no grouping)\n"
        'cp "$simulation_path/${CIF%.cif}_${SOLVENT}_net.txt" "$solvent_folder/net.txt"'
    )


def render(template: str, cfg: Config, sched: Scheduler) -> str:
    """Substitute @@TOKENS@@ for the configured values."""
    tokens = {
        "CIF": cfg.cif,
        "OCC": cfg.occ,
        "OCC_DATA": cfg.occ_data,
        "CG_EXE": cfg.cg_exe,
        "CG_KEY": cfg.cg_key,
        "VENV": cfg.venv,
        "TASK_ID": sched.task_id_var,
        "THREADS": sched.threads_var,
    }
    out = template
    for key, value in tokens.items():
        out = out.replace(f"@@{key}@@", str(value))
    return out


# ---------------------------------------------------------------------------
# Script body templates (scheduler-agnostic; tokens filled by render()).
# ---------------------------------------------------------------------------
JOB_BODY = """\
# Single OCC calculation to compute the solid-phase lattice energy.
# Run this to completion first; OCC then reuses it for every solvent.
export OCC_DATA_PATH="@@OCC_DATA@@"
occ=@@OCC@@
${occ} cg ${CIF} --radius=30 --cg-radius=3.8 --threads=${@@THREADS@@} > ${CIF%.cif}.stdout

cp ${CIF%.cif}_cg.txt STRUCTURE_FILE.txt
"""

# Common tail: lay out the solvent folder, template the input and run CG.
_SOLVENT_FOLDER_TAIL = """\
SIM_DIR=$(pwd)
cif_base=$(basename "$CIF" .cif)

# Per-solvent working folder
solvent="${SOLVENT// /+}"
simulation_path=$(realpath "$SIM_DIR")
input_path="$simulation_path/addinput.txt"

if [[ ! -f "$input_path" ]]; then
    echo "Additional input file (addinput.txt) not found."
    exit 1
fi

solvent_folder="$simulation_path/solvent_$solvent"
mkdir -p "$solvent_folder"

@@NET_BLOCK@@

ln -s "$CG_KEY" "$solvent_folder/CrystoGen.key"

# Template input.txt into the solvent folder
input_template="$simulation_path/input.txt"
input_file="$solvent_folder/input.txt"
if [[ -f "$input_template" ]]; then
    sed -e "s|FILEPATH|${solvent_folder}/|g" \\
        -e "s|SIM_NAME|${solvent}|g" \\
        "$input_template" > "$input_file"
else
    echo "Template input.txt not found in the simulation path."
    exit 1
fi

# Grow a single crystal with CrystoGen
echo "CrystoGen initiated at $(basename "$solvent_folder")"
cd "$solvent_folder"
"$CG_EXE" < "$input_path" > "${solvent_folder}/stdout.txt" 2> "${solvent_folder}/stderr.txt"
echo "Simulation Complete: $(basename "$solvent_folder")"
echo "DONE!"
"""

SCREEN_OCC_BODY = """\
@@VENV_BLOCK@@
# Environment + executables
export OCC_DATA_PATH="@@OCC_DATA@@"
OCC=@@OCC@@
@@CG_EXE_BLOCK@@

# Pick this task's solvent from solvents.txt using the array index
TASK_ID=${@@TASK_ID@@}
SOLVENT_FILE="solvents.txt"
SOLVENT=$(sed -n "${TASK_ID}p" ${SOLVENT_FILE} | awk '{$1=$1; print}')

# Compute the solvated interaction energies for this solvent
echo "Computing solvation ${CIF} in ${SOLVENT}"
${OCC} cg ${CIF} --radius=30 --cg-radius=3.8 --threads=${@@THREADS@@} --solvent="${SOLVENT}" > "${CIF%.cif}.${SOLVENT}.stdout"
echo "Done with OCC solvation computation"

""" + _SOLVENT_FOLDER_TAIL

SCREEN_NO_OCC_BODY = """\
@@VENV_BLOCK@@
# Executables (net files are assumed to already exist from a prior OCC run)
@@CG_EXE_BLOCK@@

# Pick this task's solvent from solvents.txt using the array index
TASK_ID=${@@TASK_ID@@}
SOLVENT_FILE="solvents.txt"
SOLVENT=$(sed -n "${TASK_ID}p" ${SOLVENT_FILE} | awk '{$1=$1; print}')

""" + _SOLVENT_FOLDER_TAIL

MIXED_BODY = """\
@@CG_EXE_BLOCK@@

# Pick this task's net file from file_list.txt using the array index
LIST_FILE="file_list.txt"
TASK_ID=${@@TASK_ID@@}
NET_COPY=$(sed -n "${TASK_ID}p" ${LIST_FILE} | awk '{$1=$1; print}')

SIM_DIR=$(pwd)
simulation_path=$(realpath "$SIM_DIR")
input_path="$simulation_path/addinput.txt"

if [[ ! -f "$input_path" ]]; then
    echo "Additional input file (addinput.txt) not found."
    exit 1
fi

# Numbered working folder for this mixture
job_folder="$simulation_path/$TASK_ID"
mkdir -p "$job_folder"

cp "$NET_COPY" "$job_folder/net.txt"
ln -s "$CG_KEY" "$job_folder/CrystoGen.key"

input_template="$simulation_path/input.txt"
input_file="$job_folder/input.txt"
if [[ -f "$input_template" ]]; then
    sed -e "s|FILEPATH|${job_folder}/|g" \\
        -e "s|SIM_NAME|${TASK_ID}|g" \\
        "$input_template" > "$input_file"
else
    echo "Template input.txt not found in the simulation path."
    exit 1
fi

echo "CrystoGen initiated at $(basename "$job_folder") (MIXTURE)"
cd "$job_folder"
"$CG_EXE" < "$input_path" > "${job_folder}/stdout.txt" 2> "${job_folder}/stderr.txt"
echo "Simulation Complete: $(basename "$job_folder")"
echo "DONE!"
"""

GROWTH_RATES_BODY = """\
# Usage: <submit> <script> <start> <increment>
#   start     = starting supersaturation (kcal/mol)
#   increment = step added per array task (kcal/mol)
var=$1
inc=$2

@@CG_EXE_BLOCK@@

TASK_ID=${@@TASK_ID@@}
SIM_DIR=$(pwd)
simulation_path=$(realpath "$SIM_DIR")
input_path="$simulation_path/addinput.txt"

if [[ ! -f "$input_path" ]]; then
    echo "Additional input file (addinput.txt) not found."
    exit 1
fi

job_folder="$simulation_path/$TASK_ID"
mkdir -p "$job_folder"

cp net.txt "$job_folder/net.txt"
ln -s "$CG_KEY" "$job_folder/CrystoGen.key"

input_template="$simulation_path/input.txt"
input_file="$job_folder/input.txt"
if [[ -f "$input_template" ]]; then
    var=$(echo "$var + $inc * ($TASK_ID - 1)" | bc)
    sed -e "s|FILEPATH|${job_folder}/|g" \\
        -e "s|SIM_NAME|${TASK_ID}|g" \\
        -e "s|supersat_vary|${var}|g" \\
        "$input_template" > "$input_file"
else
    echo "Template input.txt not found in the simulation path."
    exit 1
fi

echo "CrystoGen initiated at $(basename "$job_folder") (SUPERSATURATION ${var})"
cd "$job_folder"
"$CG_EXE" < "$input_path" > "${job_folder}/stdout.txt" 2> "${job_folder}/stderr.txt"
echo "Simulation Complete: $(basename "$job_folder")"
echo "DONE!"
"""

GROWTH_MOD_BODY = """\
# Usage: <submit> <script> <start> <increment>
#   start     = starting growth-modifier (poison) concentration
#   increment = step added per array task
var=$1
inc=$2

@@CG_EXE_BLOCK@@

TASK_ID=${@@TASK_ID@@}
SIM_DIR=$(pwd)
simulation_path=$(realpath "$SIM_DIR")
input_path="addinput.txt"

if [[ ! -f "$input_path" ]]; then
    echo "Additional input file (addinput.txt) not found."
    exit 1
fi

job_folder="$simulation_path/$TASK_ID"
mkdir -p "$job_folder"

cp net.txt "$job_folder/net.txt"
ln -s "$CG_KEY" "$job_folder/CrystoGen.key"

input_template="$simulation_path/input.txt"
input_file="$job_folder/input.txt"
addinput_template="$simulation_path/addinput.txt"
addinput_file="$job_folder/addinput.txt"

if [[ -f "$input_template" ]]; then
    sed -e "s|FILEPATH|${job_folder}/|g" \\
        -e "s|SIM_NAME|${TASK_ID}|g" \\
        "$input_template" > "$input_file"
else
    echo "Template input.txt not found in the simulation path."
    exit 1
fi

if [[ -f "$addinput_template" ]]; then
    var=$(echo "$var + $inc * ($TASK_ID - 1)" | bc)
    sed -e "s|poison_vary|${var}|g" \\
        "$addinput_template" > "$addinput_file"
else
    echo "Template addinput.txt not found in the simulation path."
    exit 1
fi

echo "CrystoGen initiated at $(basename "$job_folder") (POISONING ${var})"
cd "$job_folder"
"$CG_EXE" < "$input_path" > "${job_folder}/stdout.txt" 2> "${job_folder}/stderr.txt"
echo "Simulation Complete: $(basename "$job_folder")"
echo "DONE!"
"""

SIZE_STG2_BODY = """\
# Stage 2: regrow each solvent from the facets file produced in stage 1
# so that all crystal sizes are extracted.
@@CG_EXE_BLOCK@@

TASK_ID=${@@TASK_ID@@}
SOLVENT_FILE="solvents.txt"
SOLVENT=$(sed -n "${TASK_ID}p" ${SOLVENT_FILE} | awk '{$1=$1; print}')

SIM_DIR=$(pwd)
simulation_path=$(realpath "$SIM_DIR")

# The first task builds the shared COLOUR.txt template from any *all.txt facets
if [ "$TASK_ID" -eq 1 ]; then
    echo "MAKING TEMPLATE COLOUR FILE"
    colour=$(find . -type f -name "*all.txt" -not -empty -print -quit)
    lines=$(wc -l < "$colour")
    rm -f COLOUR.txt
    echo "  ${lines}" >> COLOUR.txt
    cat "$colour" >> COLOUR.txt
fi

solvent="${SOLVENT// /+}"
solvent_folder="$simulation_path/solvent_$solvent"

input_template="$simulation_path/input.txt"
input_file="$solvent_folder/input.txt"
colour_template="$simulation_path/COLOUR.txt"
colour_file="$solvent_folder/COLOUR.txt"

if [[ -f "$input_template" ]]; then
    sed -e "s|FILEPATH|${solvent_folder}/|g" \\
        -e "s|SIM_NAME|${solvent}|g" \\
        "$input_template" > "$input_file"
    echo "COPYING COLOUR FILE"
    cp "$colour_template" "$colour_file"
else
    echo "Template input.txt not found in the simulation path."
    exit 1
fi

input_path="../addinput_stg_2.txt"

echo "CrystoGen initiated at $(basename "$solvent_folder")"
cd "$solvent_folder"
"$CG_EXE" < "$input_path" > "${solvent_folder}/stdout.txt" 2> "${solvent_folder}/stderr.txt"
echo "Simulation Complete: $(basename "$solvent_folder")"
echo "DONE!"
"""

CG_SUB_BODY = """\
@@CG_EXE_BLOCK@@

# addinput.txt holds the responses CG reads from stdin
input_path="../addinput.txt"

# Run a single CrystoGen job in this folder
"$CG_EXE" < "$input_path" > "./stdout.txt" 2> "./stderr.txt"
echo "DONE!"
"""


def _assemble(cfg: Config, sched: Scheduler, *, head: str, body: str,
              cif_block: bool) -> str:
    """Glue header + (optional CIF block) + body and fill all tokens."""
    body = body.replace("@@NET_BLOCK@@", net_block(cfg))
    body = body.replace("@@VENV_BLOCK@@", venv_block(cfg))
    body = body.replace("@@CG_EXE_BLOCK@@", CG_EXE_BLOCK)
    parts = [head, ""]
    if cif_block:
        parts.append(CIF_BLOCK)
        parts.append("")
    parts.append(body)
    return render("\n".join(parts).rstrip() + "\n", cfg, sched)


# ---------------------------------------------------------------------------
# Submit-chain scripts (scheduler specific).
# ---------------------------------------------------------------------------
def submit_chain(cfg: Config, sched: Scheduler, stages: list[str]) -> str:
    """A submit.sh that fires off `stages` (script names) with dependencies."""
    lines = ["#!/usr/bin/env bash", "", CIF_BLOCK, ""]
    if sched.key == "slurm":
        prev = None
        for i, stage in enumerate(stages):
            var = f"jobid{i}" if i else "jobid"
            dep = f"--dependency=afterok:${{{prev}}} " if prev else ""
            lines.append(f"{var}=$(sbatch {dep}--parsable {stage} $CIF)")
            lines.append(f'echo "Submitted {stage} as ${{{var}}}"')
            prev = var
    else:  # sge: chain by job name with -hold_jid
        prev = None
        for i, stage in enumerate(stages):
            name = f"cg_stage{i}"
            hold = f"-hold_jid {prev} " if prev else ""
            lines.append(f"qsub -N {name} {hold}{stage} $CIF")
            lines.append(f'echo "Submitted {stage} as {name}"')
            prev = name
    return render("\n".join(lines).rstrip() + "\n", cfg, sched)


# ---------------------------------------------------------------------------
# Use-case builders. Each returns {filename: contents}.
# ---------------------------------------------------------------------------
def build_cg_sub(cfg, sched):
    head = header(cfg, sched, cores=screen_cores(cfg))
    return {"cg_folder_sub.sh": _assemble(cfg, sched, head=head,
                                          body=CG_SUB_BODY, cif_block=False)}


def build_solvent_screen_occ(cfg, sched):
    job = _assemble(cfg, sched, head=header(cfg, sched, cores=cfg.job_cores, ntasks=1),
                    body=JOB_BODY, cif_block=True)
    screen = _assemble(cfg, sched,
                       head=header(cfg, sched, cores=screen_cores(cfg), array=cfg.array),
                       body=SCREEN_OCC_BODY, cif_block=True)
    return {
        "job.sh": job,
        "screen.sh": screen,
        "submit.sh": submit_chain(cfg, sched, ["job.sh", "screen.sh"]),
    }


def build_solvent_screen_no_occ(cfg, sched):
    job = _assemble(cfg, sched, head=header(cfg, sched, cores=cfg.job_cores, ntasks=1),
                    body=JOB_BODY, cif_block=True)
    screen = _assemble(cfg, sched,
                       head=header(cfg, sched, cores=screen_cores(cfg), array=cfg.array),
                       body=SCREEN_NO_OCC_BODY, cif_block=True)
    return {
        "job.sh": job,
        "screen_no_occ.sh": screen,
        "submit.sh": submit_chain(cfg, sched, ["job.sh", "screen_no_occ.sh"]),
    }


def build_mixed_solvent_screen(cfg, sched):
    head = header(cfg, sched, cores=screen_cores(cfg), array=cfg.array)
    return {"mixed_screen.sh": _assemble(cfg, sched, head=head,
                                         body=MIXED_BODY, cif_block=False)}


def build_growth_rates(cfg, sched):
    head = header(cfg, sched, cores=screen_cores(cfg), array=cfg.array)
    return {"growth_rate_screen.sh": _assemble(cfg, sched, head=head,
                                               body=GROWTH_RATES_BODY, cif_block=False)}


def build_growth_modifiers_screen(cfg, sched):
    head = header(cfg, sched, cores=screen_cores(cfg), array=cfg.array)
    return {"growth_mod_seq.sh": _assemble(cfg, sched, head=head,
                                           body=GROWTH_MOD_BODY, cif_block=False)}


def build_size_extract_grow_twice(cfg, sched):
    job = _assemble(cfg, sched, head=header(cfg, sched, cores=cfg.job_cores, ntasks=1),
                    body=JOB_BODY, cif_block=True)
    screen = _assemble(cfg, sched,
                       head=header(cfg, sched, cores=screen_cores(cfg), array=cfg.array),
                       body=SCREEN_OCC_BODY, cif_block=True)
    stg2 = _assemble(cfg, sched,
                     head=header(cfg, sched, cores=screen_cores(cfg), array=cfg.array),
                     body=SIZE_STG2_BODY, cif_block=True)
    return {
        "job.sh": job,
        "screen.sh": screen,
        "screen_stg_2.sh": stg2,
        "submit_all.sh": submit_chain(cfg, sched,
                                      ["job.sh", "screen.sh", "screen_stg_2.sh"]),
    }


@dataclass
class UseCase:
    description: str
    builder: object
    needs: set  # config fields worth prompting for
    default_array: int = 179


USE_CASES = {
    "cg_sub": UseCase(
        "Single CrystoGen job in a folder (no OCC, no array)",
        build_cg_sub,
        {"cg_exe", "cg_key", "cores", "walltime"},
        default_array=0,
    ),
    "solvent_screen_occ": UseCase(
        "OCC lattice energy + per-solvent OCC/CG screen (job + screen + submit)",
        build_solvent_screen_occ,
        {"cif", "occ", "occ_data", "cg_exe", "cg_key", "venv",
         "cores", "job_cores", "array", "group"},
        default_array=179,
    ),
    "solvent_screen_no_occ": UseCase(
        "Per-solvent CG screen reusing existing OCC net files (job + screen + submit)",
        build_solvent_screen_no_occ,
        {"cif", "occ", "occ_data", "cg_exe", "cg_key", "venv",
         "cores", "job_cores", "array", "group"},
        default_array=179,
    ),
    "mixed_solvent_screen": UseCase(
        "Batch of CG runs from a file_list.txt of net files",
        build_mixed_solvent_screen,
        {"cg_exe", "cg_key", "cores", "array"},
        default_array=20,
    ),
    "growth_rates": UseCase(
        "Vary supersaturation across an array (start + increment args)",
        build_growth_rates,
        {"cg_exe", "cg_key", "cores", "array"},
        default_array=20,
    ),
    "growth_modifiers_screen": UseCase(
        "Vary growth-modifier (poison) concentration across an array",
        build_growth_modifiers_screen,
        {"cg_exe", "cg_key", "cores", "array"},
        default_array=20,
    ),
    "size_extract_grow_twice": UseCase(
        "OCC + screen + regrow-from-facets sequence (job + screen + stage 2)",
        build_size_extract_grow_twice,
        {"cif", "occ", "occ_data", "cg_exe", "cg_key", "venv",
         "cores", "job_cores", "array", "group"},
        default_array=179,
    ),
}


# ---------------------------------------------------------------------------
# Interactive prompting.
# ---------------------------------------------------------------------------
def ask(label: str, default):
    suffix = f" [{default}]" if default not in (None, "") else ""
    try:
        reply = input(f"{label}{suffix}: ").strip()
    except EOFError:
        reply = ""
    return reply if reply else default


def ask_choice(label: str, choices: list[str], default: str | None = None):
    print(f"\n{label}")
    for i, choice in enumerate(choices, 1):
        marker = " (default)" if choice == default else ""
        extra = ""
        if choice in USE_CASES:
            extra = f" - {USE_CASES[choice].description}"
        print(f"  {i}) {choice}{marker}{extra}")
    while True:
        reply = ask("Choose number or name", default)
        if reply in choices:
            return reply
        if reply and reply.isdigit() and 1 <= int(reply) <= len(choices):
            return choices[int(reply) - 1]
        print("  Invalid choice, try again.")


def ask_bool(label: str, default: bool) -> bool:
    reply = ask(f"{label} (y/n)", "y" if default else "n")
    return str(reply).strip().lower().startswith("y")


def interactive_fill(args) -> Config:
    """Prompt for anything not already supplied on the command line."""
    print("=== CrystoGen job-script generator (interactive) ===")
    scheduler = args.scheduler or ask_choice(
        "Scheduler:", list(SCHEDULERS), default="slurm")
    use_case = args.use_case or ask_choice(
        "Use case:", list(USE_CASES), default="solvent_screen_occ")
    mode = args.mode or ask_choice(
        "Mode:", ["parallel", "serial"], default="parallel")

    uc = USE_CASES[use_case]
    partition_default = (args.partition or
                         (DEFAULTS["partition"] if mode == "parallel" else "serial"))

    cfg = Config(
        scheduler=scheduler, mode=mode, use_case=use_case,
        output_dir=Path(args.output_dir) if args.output_dir else Path("."),
        partition=partition_default,
        array=args.array if args.array is not None else uc.default_array,
    )

    # Output dir
    cfg.output_dir = Path(ask("Output directory", str(cfg.output_dir)))

    # Prompt only for fields this use case cares about.
    prompt_map = [
        ("cif", "CIF filename", lambda v: setattr(cfg, "cif", v)),
        ("occ", "OCC executable path", lambda v: setattr(cfg, "occ", v)),
        ("occ_data", "OCC data path (OCC_DATA_PATH)", lambda v: setattr(cfg, "occ_data", v)),
        ("cg_exe", "CrystoGen executable path", lambda v: setattr(cfg, "cg_exe", v)),
        ("cg_key", "CrystoGen.key path", lambda v: setattr(cfg, "cg_key", v)),
        ("venv", "Python venv to activate (blank for none)", lambda v: setattr(cfg, "venv", v)),
        ("cores", "Cores per task (parallel array job)", lambda v: setattr(cfg, "cores", int(v))),
        ("job_cores", "Cores for the initial OCC job", lambda v: setattr(cfg, "job_cores", int(v))),
        ("array", "Array size (number of tasks)", lambda v: setattr(cfg, "array", int(v))),
    ]
    for key, label, setter in prompt_map:
        if key not in uc.needs:
            continue
        current = getattr(cfg, key)
        setter(ask(label, current))

    if "group" in uc.needs:
        cfg.group = ask_bool("Group symmetry-equivalent interactions with cg_net.py?",
                             cfg.group)

    cfg.walltime = ask("Wall-clock time limit (blank for none, "
                       "e.g. SLURM 0-12 / SGE 12:00:00)", cfg.walltime)
    return cfg


def config_from_args(args) -> Config:
    """Build a Config straight from CLI args (non-interactive path)."""
    uc = USE_CASES[args.use_case]
    mode = args.mode or "parallel"
    partition = args.partition or (DEFAULTS["partition"] if mode == "parallel" else "serial")
    cfg = Config(
        scheduler=args.scheduler,
        mode=mode,
        use_case=args.use_case,
        output_dir=Path(args.output_dir) if args.output_dir else Path("."),
        partition=partition,
        array=args.array if args.array is not None else uc.default_array,
    )
    for key in ("cif", "occ", "occ_data", "cg_exe", "cg_key", "venv",
                "cores", "job_cores", "walltime"):
        val = getattr(args, key)
        if val is not None:
            setattr(cfg, key, val)
    if args.no_venv:
        cfg.venv = ""
    if args.group is not None:
        cfg.group = args.group
    return cfg


# ---------------------------------------------------------------------------
# Output.
# ---------------------------------------------------------------------------
def write_scripts(cfg: Config) -> None:
    sched = SCHEDULERS[cfg.scheduler]
    files = USE_CASES[cfg.use_case].builder(cfg, sched)
    out_dir = cfg.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nWriting {cfg.scheduler.upper()} / {cfg.mode} / {cfg.use_case} "
          f"scripts to {out_dir}/")
    for name, contents in files.items():
        path = out_dir / name
        path.write_text(contents)
        path.chmod(0o755)
        print(f"  + {path}")
    placeholder = "YOUR_PATH" in "".join(files.values()) or cfg.cif == "MY_CIF.cif"
    if placeholder:
        print("\nNote: some placeholder values (YOUR_PATH / MY_CIF.cif) remain — "
              "edit the generated scripts or re-run supplying the real paths.")


# ---------------------------------------------------------------------------
# CLI.
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate SGE / SLURM job scripts for CrystoGen + OCC workflows.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Run with --list to see the available use cases, or --interactive "
               "to be prompted for everything.",
    )
    p.add_argument("--scheduler", choices=list(SCHEDULERS),
                   help="Target scheduler (slurm or sge).")
    p.add_argument("--mode", choices=["parallel", "serial"],
                   help="parallel (multicore) or serial. Default: parallel.")
    p.add_argument("--use-case", dest="use_case", choices=list(USE_CASES),
                   help="Which workflow to generate.")
    p.add_argument("-o", "--output-dir", dest="output_dir",
                   help="Directory to write the scripts into (default: current dir).")

    g = p.add_argument_group("paths")
    g.add_argument("--cif", help="CIF filename (default: MY_CIF.cif).")
    g.add_argument("--occ", help="OCC executable path.")
    g.add_argument("--occ-data", dest="occ_data", help="OCC_DATA_PATH value.")
    g.add_argument("--cg-exe", dest="cg_exe", help="CrystoGen executable path.")
    g.add_argument("--cg-key", dest="cg_key", help="CrystoGen.key path.")
    g.add_argument("--venv", help="Python virtualenv to activate before grouping.")
    g.add_argument("--no-venv", action="store_true",
                   help="Do not emit a venv activation line.")

    r = p.add_argument_group("resources")
    r.add_argument("--cores", type=int, help="Cores per array task (parallel mode).")
    r.add_argument("--job-cores", dest="job_cores", type=int,
                   help="Cores for the initial OCC job.")
    r.add_argument("--array", type=int, help="Array size (number of tasks).")
    r.add_argument("--partition", help="SLURM partition (default: multicore/serial by mode).")
    r.add_argument("--walltime", help="Wall-clock limit (SLURM -t / SGE h_rt).")

    grp = p.add_mutually_exclusive_group()
    grp.add_argument("--group", dest="group", action="store_true", default=None,
                     help="Group symmetry-equivalent interactions via cg_net.py (default).")
    grp.add_argument("--no-group", dest="group", action="store_false",
                     help="Copy net files directly without grouping.")

    p.add_argument("-I", "--interactive", action="store_true",
                   help="Prompt for any values not given on the command line.")
    p.add_argument("--list", action="store_true",
                   help="List available use cases and exit.")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.list:
        print("Available use cases:\n")
        for name, uc in USE_CASES.items():
            print(f"  {name}\n      {uc.description}")
        print("\nSchedulers: slurm, sge   |   Modes: parallel, serial")
        return 0

    need_interactive = args.interactive or not (args.scheduler and args.use_case)
    if need_interactive:
        cfg = interactive_fill(args)
    else:
        cfg = config_from_args(args)

    write_scripts(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
