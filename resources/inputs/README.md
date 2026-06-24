# CrystoGen input templates

Generic templates for driving a CrystoGen (formerly
CrystalGrower) simulation. Two files are read by the engine:

| File | Role |
|------|------|
| [`input.txt`](input.txt) | Fixed-layout, **structure-independent** options — required for every run. |
| [`addinput.txt`](addinput.txt) | **Structure-/mode-dependent** options — its contents depend on the choices made in `input.txt`. |

Worked, mode-specific configurations live in [`examples/`](examples) (see
[below](#examples)). Start from the templates here, substitute the
placeholders, and add/remove `addinput.txt` lines according to the rules below.

> **Tip:** for complex configurations you can generate `input.txt`,
> `addinput.txt`, and any companion files (e.g. `screw.txt`, `colour.txt`)
> interactively with the CrystoGen GUI rather than editing them by hand.

> The engine reads values **by line position**, not by keyword. Keep one value
> per line and do not reorder, insert, or delete lines except where this README
> says a line is conditional. `addinput.txt` must contain **only values** — no
> comments or blank trailing lines.

## Placeholders

`input.txt` uses `<ANGLE_BRACKET>` tokens for the things you must set per run.
Replace every one before running:

| Placeholder | Meaning |
|-------------|---------|
| `<OUTPUT_DIR>` | Directory to store simulation output. |
| `<SIM_NAME>` | File root / simulation name (used on lines 1 and 2). |
| `<STRUCTURE_FILE>` | Structure file name, e.g. `../LTA.txt`. Must be a text file. |

Everything else is a sensible default for a **normal, network-crystal,
equilibrium** run and can be edited freely.

## `input.txt` — line by line

Lines 1–27 are values; the remaining lines are a human-readable key (ignored by
the engine). The numbered descriptions in the file match this table.

| # | Value (template default) | Meaning |
|---|--------------------------|---------|
| 1 | `<OUTPUT_DIR>/<SIM_NAME>/` | Directory to store data |
| 2 | `<SIM_NAME>` | File root to store data |
| 3 | `no` | Load a checkpoint to start? (`yes`/`no`) |
| 4 | `N/A` | File path for checkpoint load |
| 5 | `yes` | Save a checkpoint at the end? (`yes`/`no`) |
| 6 | `../<STRUCTURE_FILE>.txt` | Path to structure file (must be `.txt`) |
| 7 | `normal` | Mode: `normal` / `growth_modifier` / `ordered` (also `seed_engineering`, `surface`, `screw_stress`, `nucleation`, `diffusion`) |
| 8 | `no` | Screw dislocations? (`yes`/`no`) |
| 9 | `screw.txt` | Path to screw-dislocation file |
| 10 | `1` | Checking sweeps — `1` normal, `2` to clear internal defects |
| 11 | `1` | Entropy scale factor for surface sites (0–1), `1` = no scaling |
| 12 | `net` | Framework/`tile` (e.g. zeolite) or network/`net` (molecular, ionic) |
| 13 | `yes` | Auto-calculate required memory? (`yes`/`no`) |
| 14 | `10000` | Max memory in MBytes (default 10000) |
| 15 | `25.00` | Temperature in Celsius |
| 16 | `3000000` | Number of iterations (integer) |
| 17 | `3` | Delta-mu mode (see below) |
| 18 | `no` | Excess supersaturation of any component? (always `yes` for a MOF) |
| 19 | `100.0` | Starting delta mu [kcal/mol] |
| 20 | `1` | MOVIE: number of frames |
| 21 | `1` | MOVIE: iteration at initial frame |
| 22 | `3000000` | MOVIE: iteration at final frame |
| 23 | `1000` | Delta-mu DATA: number of outputs |
| 24 | `1` | Delta-mu DATA: iteration at initial output |
| 25 | `3000000` | Delta-mu DATA: iteration at final output |
| 26 | `yes` | Visualise crystal terraces? (`yes`/`no`) |
| 27 | `colour.txt` | Path for crystal-terrace colouring |

## `addinput.txt` — conditional layout

`addinput.txt` is **assembled from `input.txt`'s choices**. The template
provided is the common case: a normal, network (`net`) crystal, no screw
dislocations, auto-calculated memory, delta-mu mode `3`, no excess
supersaturation, with terraces visualised by `find`. That produces:

| Line | Value | Source / condition |
|------|-------|--------------------|
| 1 | `yes` | Molecules grouped? (only when line 12 = `net` and not loading a checkpoint) |
| 2 | `no` | Weight multiple bonds? (`net`, new model, no checkpoint) |
| 3 | `100000` | Iterations to equilibration (delta-mu mode dependent) |
| 4 | `1500000` | Iterations at start (delta-mu mode dependent) |
| 5 | `find` | Terrace handling: `find` / `colour` (only when line 26 = `yes`) |
| 6 | `4` | Max hkl to search (only when line above = `find`) |
| 7 | `no` | Output grown crystal? (omitted for the special modes below) |

Add, remove, or reorder lines according to the active options. The blocks are
appended in this order:

1. **Seed engineering** (line 7 = `seed_engineering`): a `cut_mode`
   (`slice`/`sphere`/`rhomb`/`hemisphere`) followed by its parameters.
2. **Screw dislocations** (line 8 = `yes`): per-screw spread/axis parameters;
   requires a `screw.txt` (see line 9).
3. **Surface** (line 7 = `surface`): substrate, thickness, `h k l`, energy,
   width, and seed `x y z`.
4. **Molecular / net** (line 12 = `net`, no checkpoint): `molecules_grouped`,
   then `weighting_multiple_bonds`, then per-coordination weights if weighting.
5. **Tiles** (line 12 = `tile`): `qn_same`, per-question values, `scalings_same`,
   per-tile scalings, baseline scaling (and an ordered penalty if line 7 =
   `ordered`).
6. **Box size** (only when line 13 = `no`, i.e. memory **not** auto-calculated,
   and no checkpoint): three lines `x y z`. With auto-memory (`yes`) these are
   omitted — as in the template.
7. **Supersaturation** (all modes except `nucleation`): the iteration counts
   selected by the delta-mu mode (line 17); then, if line 18 = `yes`, the excess
   count and per-component excess values.
8. **Terraces** (line 26 = `yes`): `find` or `colour`; then either `max_hkl`
   (for `find`) or a crystal size (for `colour`).
9. **Growth modifier** (line 7 = `growth_modifier`): poison site indices, a
   terminating `0`, then attach/remain frequencies.
10. **Screw stress** (line 7 = `screw_stress`): stress radius, start, length.
11. **Nucleation** (line 7 = `nucleation`): start/end/step size, number of sims,
    and start/end/step delta mu.
12. **Output grown** appended for the standard modes (i.e. *not*
    `seed_engineering`, `nucleation`, `growth_modifier`, `surface`,
    `screw_stress`).
13. **Checkpoint box** (line 3 = `yes`): three lines `x y z`.

### Delta-mu modes (line 17)

The delta-mu mode controls which iteration counts appear in the
supersaturation block of `addinput.txt`:

| Mode | Extra `addinput.txt` lines (in order) |
|------|----------------------------------------|
| 1 | equilibration |
| 2 | *(none)* |
| 3 | equilibration, start |
| 4 | equilibration, new delta mu, start, iterations at new delta mu |
| 5 | new delta mu, start |
| 6 | equilibration, new delta mu, start, iterations at new delta mu |
| 7 | new delta mu, second delta mu, start, iterations at new delta mu |
| 8 | *(none)* |

## Companion files

Some modes reference additional files, placed alongside `input.txt`:

- `screw.txt` — screw-dislocation definitions (when line 8 = `yes`).
- `colour.txt` — terrace colouring (when line 26 = `yes`).

## Examples

The [`examples/`](examples) folder contains complete, ready-to-adapt
configurations, each a directory with its own `input.txt`, `addinput.txt`, and
any companion files:

| Example | Demonstrates |
|---------|--------------|
| `equilibrium` | Baseline equilibrium run with auto-memory and `find` terraces. |
| `equilibrium_colour` | Terraces coloured from a `colour.txt` instead of `find`. |
| `equilibrium_no_auto` | Memory not auto-calculated — explicit box size in `addinput.txt`. |
| `equilibrium_no_find` | Terraces visualised without searching hkl. |
| `low_dmu` | Low delta-mu (mode `6`) supersaturation schedule. |
| `growth_rates` | Multi-point supersaturation sweep for growth-rate data. |
| `growth_mods` | Growth-modifier (poison) runs, with and without auto-memory. |
| `screw_dislocation` | Single screw dislocation (includes `screw.txt`). |
| `two_screw_dislocation` | Two interacting screw dislocations. |
