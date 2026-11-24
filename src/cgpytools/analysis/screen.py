"""
Crystal Shape Analysis Tool

This module provides a clean, modular approach to crystal shape analysis with support for:
- Solvent screening analysis
- General shape analysis with optional energy data
- Multiple visualisation types including Zingg plots
"""

import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Literal

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import mplcursors
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
from matplotlib import cm
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch
from natsort import natsorted
from tqdm import tqdm

from cgpytools.analysis.dimers import export_key_dimers, find_dimers_folder
from cgpytools.analysis.motifs import CrystalMotifs, Pair, find_cif
from cgpytools.analysis.shape_analysis import ShapeAnalyser
from cgpytools.analysis.surfaces import process_multiple_size_files
from cgpytools.io.crystal import CrystalShape
from cgpytools.io.log import setup_logging
from cgpytools.io.net import CGNet
from cgpytools.plot.crystal_render import draw_crystal, load_shape_points
from cgpytools.plot.plot import GlobalPlotStyler, PlotTheme, build_theme

LOG = setup_logging(name="CG-ANALYSE")

# Interaction energy columns: one per net file interaction line (Int_1, Int_2, ...), or
# one per symmetry group of interactions with --group (Group_1, Group_2, ...)
ENERGY_PREFIXES = ("Int_", "Group_")


def is_energy_column(column: str) -> bool:
    """Whether ``column`` holds an interaction energy."""
    return column.startswith(ENERGY_PREFIXES)


class ShapeAnalysisConfig:
    """Configuration class for shape analysis parameters"""

    def __init__(
        self,
        save_folder: Path,
        show_plots: bool = False,
        lmax: int = 20,
        zingg_method: Literal["bounding_box", "svd"] = "svd",
        ar_limits: bool = False,
        highlight: list[str] | None = None,
        colour_by: list[str] | None = None,
        highlight_pictures: list[str] | None = None,
        input_dir: Path | None = None,
        group_energies: bool = False,
        zoom: bool = False,
    ):
        self.save_folder = save_folder
        self.show_plots = show_plots
        self.lmax = lmax
        self.zingg_method = zingg_method
        self.ar_limits = ar_limits
        # None: no extremes labelled. Empty list: the aspect ratios. Otherwise these columns.
        self.highlight = highlight
        self.colour_by = colour_by or []
        # None: no crystal pictures. Empty list: picture what --highlight/--labels pick
        # (the aspect ratio extremes if neither is set). Otherwise these columns' extremes.
        self.highlight_pictures = highlight_pictures
        # Shape file paths are stored in the results relative to this folder
        self.input_dir = input_dir
        # Plot one energy per symmetry group of interactions (Group_*) instead of per line (Int_*)
        self.group_energies = group_energies
        # Zoom the labelled/colour-by/picture Zingg plots onto the data instead of 0-1
        self.zoom = zoom

        # Create save folder if it doesn't exist
        self.save_folder.mkdir(parents=True, exist_ok=True)

    def shape_file_entry(self, path: Path) -> str:
        """Path to record in the results for a shape file.

        Relative to ``input_dir`` where possible, so the results still point
        at the shapes when the whole folder is moved, e.g. off a cluster.
        Always uses forward slashes, so the path works on any platform.
        """
        path = Path(path)
        if self.input_dir is not None:
            try:
                return path.resolve().relative_to(Path(self.input_dir).resolve()).as_posix()
            except ValueError:
                pass
        return path.resolve().as_posix()

    def resolve_shape_file(self, entry: str) -> Path:
        """Inverse of :meth:`shape_file_entry`: the shape file a results entry refers to."""
        path = Path(entry)
        if not path.is_absolute() and self.input_dir is not None:
            path = Path(self.input_dir) / path
        return path


class FileDiscovery:
    """Handles discovery and organization of input files"""

    @staticmethod
    def find_files(
        folderpath: Path,
        include_patterns: dict[str, list[str]] | None = None,
        exclude_patterns: dict[str, list[str]] | None = None,
    ) -> dict[str, list[Path]]:
        """Discover all relevant files in the input directory"""

        # defaults
        if include_patterns is None:
            include_patterns = {
                "xyz": ["*Aspects.XYZ", "*CGvisualiser.XYZ"],
                "wulff": ["*.ply"],
                "cda": ["*simulation_parameters.txt"],
                "occ": ["*.*.stdout"],
                "size": ["*size.csv"],
            }
        if exclude_patterns is None:
            exclude_patterns = {"xyz": ["*ungrown*.XYZ"]}

        file_types: dict[str, list[Path]] = {}

        for category, patterns in include_patterns.items():
            found: list[Path] = []
            for pat in patterns:
                found.extend(folderpath.rglob(pat))

            # Apply category-specific exclusions if present
            excludes = exclude_patterns.get(category, [])
            filtered = [f for f in found if not any(f.match(excl) for excl in excludes)]

            # Natural sort
            file_types[category] = natsorted(filtered, key=lambda x: str(x))

        message = (
            f"Found:\n"
            f" {len(file_types['xyz'])} XYZs\n"
            f" {len(file_types['wulff'])} Wulff shapes\n"
            f" {len(file_types['cda'])} CDA files\n"
            f" {len(file_types['occ'])} OCC outputs\n"
            f" {len(file_types['size'])} size files\n"
        )
        LOG.info(message)

        return file_types


class EnergyDataLoader:
    """Handles loading of energy data from various sources"""

    @staticmethod
    def parse_shape_name(value: str | float) -> str | int:
        if isinstance(value, (int, float)):
            # Collapse whole floats to int
            return str(int(value) if float(value).is_integer() else str(value))
        try:
            num = float(value)
            return str(int(num) if num.is_integer() else value)
        except (ValueError, TypeError):
            return str(value)

    @staticmethod
    def load_from_csv(
        csv_path: Path, shape_name_column: str = "shape_name"
    ) -> dict[str, list[tuple[str, float]]]:
        """Load energy data from CSV file for general shape analysis"""
        try:
            df = pd.read_csv(csv_path)
            energy_data = {}

            # Get energy columns (Int_* or x*)
            energy_cols = [col for col in df.columns if col.startswith(("Int_", "x", "X"))]

            for _, row in df.iterrows():
                shape_name = EnergyDataLoader.parse_shape_name(row[shape_name_column])
                energies = []
                for col in energy_cols:
                    try:
                        energy = float(row[col])
                        if not np.isnan(energy):
                            energies.append((col, energy))
                    except (ValueError, TypeError):
                        continue
                energy_data[shape_name] = energies

            LOG.info("Loaded energy data for %s shapes from %s", len(energy_data), csv_path)
            return energy_data

        except Exception:
            LOG.exception("Failed to load energy data from %s", csv_path)
            return {}

    @staticmethod
    def load_from_csvs(
        csv_paths: list[Path], shape_name_column: str = "shape_name"
    ) -> dict[str, list[tuple[str, float]]]:
        """Load and merge energy data from multiple CSV files"""
        merged: dict[str, list[tuple[str, float]]] = {}
        for p in csv_paths or []:
            if p and p.exists():
                merged.update(EnergyDataLoader.load_from_csv(p, shape_name_column))
        return merged

    @staticmethod
    def net_file_for(shape: Path) -> Path:
        """The net.txt file holding the energies a shape was grown with."""
        shape = Path(shape)
        if shape.suffix == ".stdout":
            return shape.parent / f"solvent_{shape.parent.name.split('_')[-1]}" / "net.txt"
        return shape.parent / "net.txt"

    @staticmethod
    def find_net_file(df: pd.DataFrame, config: "ShapeAnalysisConfig") -> Path | None:
        """A net file behind the results in ``df``, to read the energy groups from.

        Tried first beside the shape files the results record, then anywhere
        under the input folder. Every net file in a screen describes the same
        crystal, so any one of them gives the groups.
        """
        if "shape_file" in df.columns:
            for entry in df["shape_file"].dropna().unique():
                netfile = EnergyDataLoader.net_file_for(config.resolve_shape_file(entry))
                if netfile.exists():
                    return netfile
        if config.input_dir is not None and Path(config.input_dir).is_dir():
            return next(Path(config.input_dir).rglob("net.txt"), None)
        return None

    @staticmethod
    def read_energy_groups(netfile: Path) -> dict[str, list[dict]] | None:
        """The interaction lines of a net file, grouped by their shared ID.

        Lines related by symmetry share an ID, and so one energy. Each group is
        named ``Group_<ID>`` (``Group_<label>_<ID>`` if IDs repeat across
        molecules) and lists its lines in file order as dicts with ``column``
        (the ``Int_N`` column of that line), ``dimer`` (e.g. ``dimer_6``),
        ``r`` (the centroid distance) and ``pair`` (the line's two molecules as
        a :class:`~cgpytools.analysis.motifs.Pair`, or None). Returns None if
        the file can't be read.
        """
        try:
            net = CGNet(netfile)
            net.parse()
            interactions = net.energy_interactions()
            dimers = net.energy_dimers
        except Exception:
            LOG.exception("--group: failed to read the energy groups from %s", netfile)
            return None

        labels_per_id = defaultdict(set)
        for label, interaction in interactions:
            labels_per_id[interaction.id].add(label)
        with_label = any(len(labels) > 1 for labels in labels_per_id.values())

        groups: dict[str, list[dict]] = {}
        for i, ((label, interaction), dimer) in enumerate(zip(interactions, dimers), 1):
            name = f"Group_{label}_{interaction.id}" if with_label else f"Group_{interaction.id}"
            groups.setdefault(name, []).append(
                {
                    "column": f"Int_{i}",
                    "dimer": dimer,
                    "r": interaction.r,
                    "pair": Pair.parse(interaction.molecule_info, interaction.r),
                }
            )
        return groups

    @staticmethod
    def group_energy_columns(
        df: pd.DataFrame, groups: dict[str, list[dict]], netfile: Path
    ) -> pd.DataFrame:
        """Replace the per-line ``Int_*`` columns with one column per group of
        :meth:`read_energy_groups`.

        ``Int_N`` is the energy of the Nth interaction line of the net file.
        Each group column holds the mean of its lines, which are equal unless
        the net file was edited by hand. ``df`` is returned unchanged if the
        net file does not match its energy columns.
        """
        int_cols = natsorted(col for col in df.columns if re.fullmatch(r"Int_\d+", col))
        n_lines = sum(len(lines) for lines in groups.values())
        if int_cols != [f"Int_{i}" for i in range(1, n_lines + 1)]:
            LOG.warning(
                "--group: %s has %d interaction lines but the results have %d energy columns "
                "(%s ... %s); plotting ungrouped energies",
                netfile,
                n_lines,
                len(int_cols),
                int_cols[0],
                int_cols[-1],
            )
            return df
        groups = {name: [line["column"] for line in lines] for name, lines in groups.items()}

        grouped = {}
        for name, members in groups.items():
            values = df[members].apply(pd.to_numeric, errors="coerce")
            spread = (values.max(axis=1) - values.min(axis=1)).max()
            if spread > 1e-6:
                LOG.warning(
                    "--group: %s energies differ by up to %.4g; using their mean",
                    ", ".join(members),
                    spread,
                )
            grouped[name] = values.mean(axis=1)

        LOG.info("Energy groups from %s:", netfile)
        for name, members in groups.items():
            LOG.info("  %-12s %s", name, ", ".join(members))

        # Grouped columns take the place of the first Int_ column, keeping the column order
        position = df.columns.get_loc(int_cols[0])
        result = df.drop(columns=int_cols)
        for offset, (name, values) in enumerate(grouped.items()):
            result.insert(position + offset, name, values)
        return result

    @staticmethod
    def load_from_net_files(shapes: list[Path]) -> dict[str, list[float]]:
        """Load energy data from net.txt files"""
        energy_data = {}

        for shape in shapes:
            shape = Path(shape)
            netfile = EnergyDataLoader.net_file_for(shape)

            if netfile.exists():
                try:
                    net = CGNet(netfile)
                    net.parse()
                    energies = net.unique_energies_arr.flatten()
                    energies = energies[~np.isnan(energies)]  # Remove NaN values
                    energy_data[shape.stem] = energies.tolist()
                    LOG.debug("Loaded %s energies from %s", len(energies), netfile)
                except Exception:
                    LOG.exception("Failed to load energies from %s", netfile)

        return energy_data


class SolventDataLoader:
    """Handles loading and processing of solvent-related data"""

    @staticmethod
    def load_solvent_properties(solvent_json: Path) -> dict:
        """Load solvent properties from JSON file"""

        try:
            if not Path(solvent_json).is_file or Path(solvent_json).stat().st_size == 0:
                LOG.error("File is empty/not found: %s", solvent_json)
                return FileNotFoundError
            with open(solvent_json, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            LOG.exception("Failed to load solvent properties from %s", solvent_json)
            return {}

    @staticmethod
    def extract_solvent_name(file_path: Path) -> str | None:
        """Extract solvent name from file path or filename"""
        name_split = str(file_path.parent).rsplit("_", maxsplit=1)
        solvent = name_split[-1] if len(name_split) > 1 else None

        if solvent:
            # Clean up solvent names
            solvent = solvent.replace("+", " ")
            solvent = solvent.replace("E-Z", "E/Z")
            solvent = solvent.replace("cis-trans", "cis/trans")

        return solvent

    @staticmethod
    def load_occ_solubility(occ_outputs: list[Path]) -> dict[str, list[float]]:
        """Load solubility data from OCC output files"""
        occ_info = defaultdict(list)

        for output in occ_outputs:
            output = Path(output)
            solvent = output.name.split(".")[-2]

            try:
                with open(output, "r", encoding="utf-8") as f:
                    lines = f.readlines()

                solubility = None
                for line in lines:
                    if line.startswith("solubility (g/L)"):
                        solubility = float(line.split()[-1])
                        break

                if solubility is not None:
                    LOG.info("Solubility in %s: %s", solvent, solubility)
                    occ_info[solvent].append(solubility)
                else:
                    LOG.warning("Solubility not found in %s", output)

            except Exception:
                LOG.exception("Failed to process OCC output %s", output)

        return dict(occ_info)


class CrystalShapeAnalyser:
    """Main class for analyzing crystal shapes"""

    def __init__(self, config: ShapeAnalysisConfig):
        self.config = config
        self.shape_analyser = ShapeAnalyser(zingg_method=self.config.zingg_method)

    def analyse_general_shapes(
        self,
        shape_files: list[Path],
        energy_csv: Path | None = None,
    ) -> pd.DataFrame:
        """Analyse general shapes with optional energy data from CSV"""
        shape_info = defaultdict(list)
        energy_data = {}

        if energy_csv and energy_csv.exists():
            energy_data = EnergyDataLoader.load_from_csv(energy_csv)

        for shape_file in tqdm(shape_files, desc="Analysing general shapes"):
            shape_file = Path(shape_file)
            shape_name = shape_file.stem.replace("_CGvisualiser", "").replace("_CGAspects", "")

            try:
                crystal = CrystalShape.from_file(shape_file)
                self.shape_analyser.analyse_crystal(crystal, frame_idx=-1)
                metrics = self.shape_analyser.get_frame_metrics(-1)

                shape_info["shape_name"].append(shape_name)
                shape_info["shape_file"].append(self.config.shape_file_entry(shape_file))
                shape_info["ar1"].append(metrics.aspect1)
                shape_info["ar2"].append(metrics.aspect2)
                shape_info["sa"].append(metrics.surface_area)
                shape_info["vol"].append(metrics.volume)
                shape_info["sa_vol"].append(metrics.surface_area_to_volume_ratio)

                if shape_name in energy_data:
                    for i, (col, energy) in enumerate(energy_data[shape_name], 1):
                        shape_info[col].append(energy)

                LOG.debug("Analysed shape: %s", shape_name)

            except Exception:
                LOG.exception("Failed to analyse shape %s", shape_file)
                continue
        LOG.info("\n")

        return self._finalise_dataframe(shape_info, "general_shapes")

    def analyse_movies(
        self,
        movie_files: list[Path],
        energy_csv: Path | None = None,
    ) -> pd.DataFrame:
        """Analyse shape movies (multi-frame shapes) with optional energy data from CSV."""
        movie_info = defaultdict(list)
        energy_data = {}

        if energy_csv and energy_csv.exists():
            energy_data = EnergyDataLoader.load_from_csv(energy_csv)

        for movie_file in tqdm(movie_files, desc="Analysing shape movies"):
            movie_file = Path(movie_file)
            movie_name = movie_file.stem.replace("_CGvisualiser", "").replace("_CGAspects", "")

            try:
                # Load the multi-frame crystal shape
                crystal = CrystalShape.from_file(movie_file)

                # Analyse all frames in the movie
                self.shape_analyser.analyse_crystal(crystal, frame_idx=None)

                movie_entry = self.config.shape_file_entry(movie_file)
                for frame_idx, metrics in self.shape_analyser.get_all_frame_metrics().items():
                    movie_info["shape_name"].append(movie_name)
                    movie_info["shape_file"].append(movie_entry)
                    movie_info["frame"].append(frame_idx)
                    movie_info["ar1"].append(metrics.aspect1)
                    movie_info["ar2"].append(metrics.aspect2)
                    movie_info["sa"].append(metrics.surface_area)
                    movie_info["vol"].append(metrics.volume)
                    movie_info["sa_vol"].append(metrics.surface_area_to_volume_ratio)

                    # Add energy values if available (same mapping as general shapes)
                    if movie_name in energy_data:
                        for col, energy in energy_data[movie_name]:
                            movie_info[col].append(energy)

                LOG.debug("Analysed movie: %s (%s frames)", movie_name, len(crystal.frames))

            except Exception:
                LOG.exception("Failed to analyse movie %s", movie_file)
                continue
        LOG.info("\n")

        return self._finalise_dataframe(movie_info, "movies")

    def analyse_solvent_shapes(
        self,
        shapes: list[Path],
        solvent_json: Path,
        occ_outputs: list[Path] | None = None,
        get_energies: bool = False,
    ) -> pd.DataFrame:
        """Analyse shapes from solvent screening"""
        shape_info = defaultdict(list)

        sol_dict = SolventDataLoader.load_solvent_properties(solvent_json)
        occ_info = SolventDataLoader.load_occ_solubility(occ_outputs) if occ_outputs else {}
        energy_data = EnergyDataLoader.load_from_net_files(shapes) if get_energies else {}

        for shape_file in tqdm(shapes, desc="Analysing solvent shapes"):
            shape_file = Path(shape_file)
            solvent = self._get_solvent_from_shape(shape_file)
            if not solvent or solvent not in sol_dict:
                LOG.error("Couldn't find solvent for %s", shape_file)
                continue
            try:
                crystal = CrystalShape.from_file(shape_file)
                self.shape_analyser.analyse_crystal(
                    crystal,
                    frame_idx=-1,
                )
                metrics = self.shape_analyser.get_frame_metrics(-1)
                shape_info["solvent"].append(solvent)
                shape_info["shape_file"].append(self.config.shape_file_entry(shape_file))
                shape_info["ar1"].append(metrics.aspect1)
                shape_info["ar2"].append(metrics.aspect2)
                shape_info["sa"].append(metrics.surface_area)
                shape_info["vol"].append(metrics.volume)
                shape_info["sa_vol"].append(metrics.surface_area_to_volume_ratio)
                params = sol_dict[solvent]
                shape_info["n"].append(params[0])
                shape_info["acidity"].append(params[1])
                shape_info["basicity"].append(params[2])
                shape_info["gamma"].append(params[3])
                shape_info["dielectric"].append(params[4])
                shape_info["aromatic"].append(params[5])
                shape_info["halogen"].append(params[6])
                if solvent in occ_info:
                    shape_info["solubility"].append(np.log10(occ_info[solvent][0]))
                shape_stem = shape_file.stem
                if shape_stem in energy_data:
                    try:
                        for i, (col, energy) in enumerate(energy_data[shape_stem], 1):
                            shape_info[col].append(energy)
                    except TypeError:
                        for i, energy in enumerate(energy_data[shape_stem], 1):
                            shape_info[f"Int_{i}"].append(energy)
                tqdm.write(f"Analysed shape from solvent: {solvent}")
            except Exception:
                LOG.exception("Failed to analyse shape %s", shape_file)
                continue

        LOG.info("\n")

        return self._finalise_dataframe(shape_info, "solvent_shapes")

    def analyse_cda_shapes(
        self, cda_files: list[Path], directions: list[str], get_energies: bool = False
    ) -> pd.DataFrame:
        """Analyse CDA simulation results"""
        ar_dict = defaultdict(list)

        for cda_file in tqdm(cda_files, desc="Analyzing CDA files"):
            solvent = SolventDataLoader.extract_solvent_name(cda_file)
            if not solvent:
                LOG.warning("Could not extract solvent name from %s", cda_file)
                continue

            try:
                with open(cda_file, "r", encoding="utf-8") as f:
                    lines = f.readlines()

                frame_start = None
                for i, line in enumerate(lines):
                    if line.startswith("Size of crystal at frame output"):
                        frame_start = i + 1
                        break

                if frame_start is None:
                    LOG.warning("Frame information not found in %s", cda_file)
                    continue

                len_info_lines = lines[frame_start:]
                for direction in directions:
                    for len_line in len_info_lines:
                        if len_line.startswith(direction):
                            ar_dict[direction].append(float(len_line.split(" ")[-2]))
                            break

                if solvent not in ar_dict["solvent"]:
                    ar_dict["solvent"].append(solvent)

                if get_energies:
                    netfile = cda_file.parent / "net.txt"
                    if netfile.exists():
                        try:
                            net = CGNet(netfile)
                            net.parse()
                            energies = net.unique_energies_arr.flatten()
                            energies = energies[~np.isnan(energies)]
                            for i, energy in enumerate(energies, 1):
                                ar_dict[f"Int_{i}"].append(energy)
                        except Exception:
                            LOG.exception("Failed to load energies from %s", netfile)

            except Exception:
                LOG.exception("Failed to process CDA file %s", cda_file)
                continue

        df = pd.DataFrame.from_dict(ar_dict)
        for i in range(len(directions) - 1):
            df[f"AspectRatio_{directions[i]}/{directions[i + 1]}"] = (
                df[directions[i]] / df[directions[i + 1]]
            )

        df.to_csv(self.config.save_folder / "cda_analysis.csv", index=False)
        LOG.info("\n")
        return df

    def analyse_size_wulff_shapes(
        self,
        generation_results: dict[str, dict],
        energy_csvs: list[Path] | None = None,
        relative: bool = True,
    ) -> pd.DataFrame:
        """
        Analyze generated Wulff shapes from size files using the existing ShapeAnalyser.

        Args:
            generation_results: Results from generation step
            energy_csvs: Optional list of CSV paths with energy data (one per subfolder group)
            relative: whether surface area and volume are zeroed at initial frame

        Returns:
            DataFrame with shape analysis results. If size files are organised in
            subdirectories, 'subfolder' (parent.parent dir name) and 'group'
            (parent.parent.parent dir name) columns are added when they vary across files.
        """

        shape_info = defaultdict(list)

        # Pre-compute path-derived labels to decide which columns are meaningful
        all_subfolders = {Path(sf).parent.parent.name for sf in generation_results}
        all_groups = {Path(sf).parent.parent.parent.name for sf in generation_results}
        add_subfolder = len(all_subfolders) >= 2
        add_group = len(all_groups) >= 2

        # Load energy data if provided
        energy_data = EnergyDataLoader.load_from_csvs(energy_csvs)

        for size_file, size_rows in tqdm(
            generation_results.items(), desc="Analyzing Size Wulff shapes"
        ):
            size_name = Path(size_file).stem.replace("_size", "")
            subfolder = Path(size_file).parent.parent.name if add_subfolder else None
            group = Path(size_file).parent.parent.parent.name if add_group else None
            sa_min = 0
            vol_min = 0
            for i, r in enumerate(size_rows):
                timestep = r["timestep"]
                shape_file = Path(r["output_file"])

                try:
                    # Analyze the shape using existing infrastructure
                    crystal = CrystalShape.from_file(shape_file)
                    # Use the existing shape analyser
                    self.shape_analyser.analyse_crystal(
                        crystal,
                        frame_idx=-1,
                    )
                    metrics = self.shape_analyser.get_frame_metrics(-1)

                    # Store basic shape information
                    shape_info["shape_name"].append(size_name)
                    shape_info["shape_file"].append(self.config.shape_file_entry(shape_file))
                    shape_info["timestep"].append(timestep)
                    shape_info["time"].append(r.get("time_value", timestep))
                    if i == 0:
                        sa_min = metrics.surface_area
                        vol_min = metrics.volume

                    shape_info["ar1"].append(metrics.aspect1)
                    shape_info["ar2"].append(metrics.aspect2)
                    shape_info["sa"].append(
                        metrics.surface_area - sa_min if relative else metrics.surface_area
                    )
                    shape_info["vol"].append(
                        metrics.volume - vol_min if relative else metrics.volume
                    )
                    shape_info["sa_vol"].append(metrics.surface_area_to_volume_ratio)

                    if group is not None:
                        shape_info["group"].append(group)
                    if subfolder is not None:
                        shape_info["subfolder"].append(subfolder)

                    # Add energy data if available
                    if size_name in energy_data:
                        for col, energy in energy_data[size_name]:
                            shape_info[col].append(energy)

                    LOG.debug("Analyzed Wulff shape: %s", shape_file.name)

                except Exception:
                    LOG.exception("Failed to analyze shape %s at row %d", size_file, i)
                    continue

        if not shape_info:
            LOG.warning("No shapes were successfully analyzed")
            return pd.DataFrame()

        # Convert to DataFrame and return
        df = pd.DataFrame(shape_info)
        LOG.info("Successfully analyzed %s Wulff shapes", len(df))
        df.to_csv(self.config.save_folder / "size_analysis.csv", index=False)

        return df

    def _get_solvent_from_shape(self, shape_file: Path) -> str | None:
        """Extract solvent information from shape file"""
        # Try to get from parent directory
        solvent = SolventDataLoader.extract_solvent_name(shape_file)

        # If not found, try to read from SHAPE file
        if not solvent and shape_file.name.startswith("SHAPE"):
            try:
                with open(shape_file, "r", encoding="utf-8") as f:
                    lines = f.readlines()

                if len(lines) > 1:
                    sol_line = lines[1]
                    if sol_line.startswith("Solvent"):
                        solvent = (
                            sol_line.split(":")[-1]
                            .strip()
                            .replace("\n", "")
                            .replace("E-Z", "E/Z")
                            .replace("cis-trans", "cis/trans")
                        )
            except Exception:
                LOG.exception("Failed to read solvent from %s", shape_file)

        return solvent

    def _finalise_dataframe(self, shape_info: dict, name: str) -> pd.DataFrame:
        """Finalise the shape information dictionary into a DataFrame"""
        # Find the maximum list length
        max_length = max(len(lst) for lst in shape_info.values() if lst)

        # Pad shorter lists with zeros and warn
        for key, value in shape_info.items():
            if len(value) < max_length:
                original_length = len(value)
                shape_info[key] = value + [0] * (max_length - len(value))
                LOG.warning(
                    "%s values were padded with %s zeros",
                    key,
                    max_length - original_length,
                )

        df = pd.DataFrame(shape_info)
        df.to_csv(self.config.save_folder / f"{name}.csv", index=False)

        LOG.info(df)
        LOG.info("\n")

        return df


class ShapePlots:
    """Handles all visualisation tasks with consistent global styling"""

    def __init__(
        self,
        config: ShapeAnalysisConfig,
        plot_style: Literal["modern", "classic", "minimal", "dark", "publication"] = "modern",
        custom_theme: PlotTheme = None,
    ):
        self.config = config

        # Initialize global styling
        self.styler = GlobalPlotStyler(theme=custom_theme, style=plot_style)

        # Apply column renaming function for consistent naming
        self._clean_columns = lambda df: df.rename(
            columns=lambda x: x.replace("x_", "").replace("X_", "")
        )

    @staticmethod
    def get_subplot_layout(n: int) -> tuple[int, int]:
        """Calculate optimal subplot layout"""
        n_cols = math.ceil(math.sqrt(n))
        n_rows = math.ceil(n / n_cols)
        return n_rows, n_cols

    def _apply_ar_limits(self, ax, x_is_ar=True, y_is_ar=True):
        """Apply aspect ratio limits if configured"""
        if self.config.ar_limits:
            if x_is_ar:
                ax.set_xlim(0, 1)
            if y_is_ar:
                ax.set_ylim(0, 1)

    @staticmethod
    def _stack_label_slots(
        fracs: list[float], min_gap: float, lo: float = 0.04, hi: float = 0.96
    ) -> list[float]:
        """Push overlapping label positions apart, keeping them inside [lo, hi].

        ``fracs`` are the preferred y positions in axes fraction, one per label.
        Returns positions in the same order, each at least ``min_gap`` from its
        neighbour, so labels for points in the same region do not collide.
        """
        order = sorted(range(len(fracs)), key=lambda i: fracs[i])

        placed: list[float] = []
        for i in order:
            target = min(max(fracs[i], lo), hi)
            if placed:
                target = max(target, placed[-1] + min_gap)
            placed.append(target)

        # If the stack ran off the top, slide it down and re-separate downwards
        overflow = placed[-1] - hi
        if overflow > 0:
            placed = [frac - overflow for frac in placed]
            for i in range(len(placed) - 2, -1, -1):
                placed[i] = min(placed[i], placed[i + 1] - min_gap)

        slots = [0.0] * len(fracs)
        for frac, i in zip(placed, order):
            slots[i] = frac
        return slots

    def _place_leader_labels(
        self,
        ax,
        points: list[tuple[float, float]],
        texts: list[str],
        box_color: str,
        arrow_color: str,
    ):
        """Label ``points`` from a column of boxes just outside the right of the axes.

        Keeping the boxes outside the plot means they never cover data points or
        other leader lines. Boxes are ordered top to bottom like their points
        and de-cluttered by :meth:`_stack_label_slots`; the leader lines are
        drawn beneath highlighted markers (zorder above 4) so those stay visible.
        Call after the axes limits are final.
        """
        y0, y1 = ax.get_ylim()
        y_span = (y1 - y0) or 1.0

        theme = self.styler.theme
        # A label box is roughly 2.6x the font size tall, in axes fractions
        fig_height_pts = ax.figure.get_size_inches()[1] * 72
        min_gap = max(2.6 * theme.font_size_tick / fig_height_pts, 0.035)

        items = [
            (point_x, point_y, text)
            for (point_x, point_y), text in zip(points, texts)
            if np.isfinite(point_x) and np.isfinite(point_y)
        ]
        if not items:
            return

        slots = self._stack_label_slots(
            [(point_y - y0) / y_span for _, point_y, _ in items], min_gap, lo=0.02, hi=0.98
        )

        for (point_x, point_y, text), frac in zip(items, slots):
            ax.annotate(
                text,
                xy=(point_x, point_y),
                xycoords="data",
                xytext=(1.03, frac),
                textcoords="axes fraction",
                ha="left",
                va="center",
                bbox={
                    "boxstyle": "round,pad=0.35",
                    "fc": box_color,
                    "alpha": 0.95,
                    "edgecolor": arrow_color,
                    "linewidth": 0.8,
                },
                arrowprops={
                    "arrowstyle": "-",
                    "color": arrow_color,
                    "alpha": 0.8,
                    "linewidth": 1.0,
                    "shrinkA": 0,
                    "shrinkB": 5,
                },
                fontsize=theme.font_size_tick,
                fontweight="medium",
                annotation_clip=False,
                zorder=4,
            )

    @staticmethod
    def _find_extremes(columns: dict[str, pd.Series], warn: bool = True) -> dict[int, list[str]]:
        """Find the rows holding the minimum and maximum of each column.

        ``columns`` maps a display name (e.g. ``"S:M"``) to its values. Returns
        a mapping of row index to tags such as ``["max S:M", "min M:L"]``, so a
        row that is extreme in more than one column gets a single label.
        Ties resolve to the first row; non-numeric and NaN values are ignored.
        """
        extremes: dict[int, list[str]] = defaultdict(list)
        for display_name, values in columns.items():
            numeric = pd.to_numeric(values, errors="coerce").dropna()
            if numeric.empty:
                if warn:
                    LOG.warning(
                        "No numeric values in %s, cannot highlight its extremes", display_name
                    )
                continue
            extremes[numeric.idxmin()].append(f"min {display_name}")
            extremes[numeric.idxmax()].append(f"max {display_name}")
        return dict(extremes)

    def _zoom_to_data(self, ax, x: pd.Series, y: pd.Series, pad: float = 0.04):
        """Tighten Zingg axes onto the data if --zoom was given.

        Zingg plots keep their 0-1 axes by default, so the shape regions stay
        in proportion. ``apply_zingg_style`` always sets 0-1 limits, so this
        must be called after it to have any effect.
        """
        if not self.config.zoom or self.config.ar_limits:
            return

        finite = np.isfinite(x) & np.isfinite(y)
        if not finite.any():
            return

        for setter, values in ((ax.set_xlim, x[finite]), (ax.set_ylim, y[finite])):
            lo, hi = float(values.min()), float(values.max())
            margin = max((hi - lo) * pad, pad / 2)
            setter(max(0.0, lo - margin), min(1.0, hi + margin))

    def create_zingg_plot(self, df: pd.DataFrame, name: str = "", interactive: bool = True):
        """Create basic Zingg classification plot with consistent styling"""
        # Clean column names
        df_clean = self._clean_columns(df)

        fig, ax = self.styler.create_figure(figsize=(10, 8))

        # Determine which columns to use for aspect ratios
        if name == "cda":
            ar_cols = [col for col in df_clean.columns if col.startswith("AspectRatio")]
            if len(ar_cols) < 2:
                LOG.error("Not enough aspect ratio columns for CDA plot")
                return
            x, y = df_clean[ar_cols[0]], df_clean[ar_cols[1]]
            xlabel, ylabel = ar_cols[0], ar_cols[1]
        else:
            x, y = df_clean["ar1"], df_clean["ar2"]
            xlabel, ylabel = "S:M", "M:L"
            self._apply_ar_limits(ax, x_is_ar=True, y_is_ar=True)

        LOG.info("Plotting Basic Zingg Plot (%s)", name)

        # Create scatter plot with styled colors
        colors = self.styler.get_color_palette(1)[0]
        scatter = ax.scatter(
            x,
            y,
            alpha=self.styler.theme.alpha_scatter,
            s=self.styler.theme.marker_size,
            color=colors,
        )

        # Apply Zingg styling (includes classification lines and regions)
        title = (
            f"Zingg Classification Plot - {name.upper()}" if name else "Zingg Classification Plot"
        )
        self.styler.apply_zingg_style(ax, title=title, show_legend=True)

        # Override xlabel/ylabel for CDA plots
        if name == "cda":
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

        # Add interactive tooltips
        if hasattr(df_clean, "solvent"):
            tooltip_data = df_clean["solvent"]
        elif hasattr(df_clean, "shape_name"):
            tooltip_data = df_clean["shape_name"]
        else:
            tooltip_data = df_clean.index

        cursors = mplcursors.cursor(scatter, hover=True)
        cursors.connect(
            "add", lambda sel: sel.annotation.set_text(str(tooltip_data.iloc[sel.target.index]))
        )

        plt.tight_layout()

        # Save with consistent settings
        save_path = self.config.save_folder / f"zingg_{name}.png"
        self.styler.save_figure(fig, save_path)

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

        # Create interactive version
        if interactive:
            fig_interactive = px.scatter(
                x=x,
                y=y,
                hover_name=tooltip_data,
                title=f"Interactive Zingg Plot - {name.upper()}",
                color_discrete_sequence=self.styler.theme.categorical_palette,
            )
            if self.config.ar_limits:
                fig_interactive.update_xaxes(range=[0, 1])
                fig_interactive.update_yaxes(range=[0, 1])

            fig_interactive.write_html(self.config.save_folder / f"zingg_{name}_interactive.html")

    def plot_var_heatmaps(
        self,
        df: pd.DataFrame,
        name: str = "",
        mode="parameter",
        z_vars=None,
        timecol: Literal["timestep", "time", "frame"] = "timestep",
        smooth_window: int = 1,
        log_scale: bool = False,
    ):
        """
        Heatmap subplots with all x_* variables vs specified z variables

        Args:
            df: DataFrame with data
            name: Name for saving
            mode: "parameter" or "energy"
            z_vars: List of z-axis variables, or None to use default ["sa_vol"]
        """
        # Clean column names first
        df_clean = self._clean_columns(df)

        # Default z variables if none specified
        if z_vars is None:
            z_vars = ["sa_vol"]

        # Find x variables (parameter columns)
        if mode == "energy":
            # Energy columns
            x_cols = [col for col in df_clean.columns if is_energy_column(col)]
            if not x_cols:
                LOG.warning("No energy columns found")
                return
        elif mode == "parameter":
            # Exclude common/energy columns, use the rest
            common_cols = [
                "solvent",
                "shape_name",
                "shape_file",
                "name",
                "ar1",
                "ar2",
                "sa",
                "vol",
                "sa_vol",
                "timestep",
                "time",
                "frame",
                "subfolder",
                "group",
            ]
            int_cols = [col for col in df_clean.columns if is_energy_column(col)]
            exclude = set(common_cols + int_cols)
            x_cols = [col for col in df_clean.columns if col not in exclude]
        else:
            LOG.error("Unknown mode '%s'. Use 'energy' or 'parameter'.", mode)
            return

        if not x_cols:
            LOG.warning("No x columns found for mode=%s", mode)
            return

        # Check that all z_vars exist in the dataframe
        missing_z_vars = [z for z in z_vars if z not in df_clean.columns]
        if missing_z_vars:
            LOG.error("Missing z variables in dataframe: %s", missing_z_vars)
            return

        # Create separate figures for each z variable
        has_subfolders = "subfolder" in df_clean.columns
        for z_var in z_vars:
            LOG.info("Plotting Heatmap Plot for %s Variables vs %s (%s)", mode.title(), z_var, name)
            if has_subfolders:
                unique_subfolders = sorted(df_clean["subfolder"].dropna().unique())
                for x_var in x_cols:
                    self._create_subfolder_heatmap_figure(
                        df_clean,
                        x_var,
                        unique_subfolders,
                        z_var,
                        name,
                        timecol,
                        smooth_window=smooth_window,
                        log_scale=log_scale,
                    )
            else:
                self._create_heatmap_figure(
                    df_clean,
                    x_cols,
                    z_var,
                    name,
                    mode,
                    timecol,
                    smooth_window=smooth_window,
                    log_scale=log_scale,
                )

    def _create_subfolder_heatmap_figure(
        self,
        df_clean,
        x_var,
        unique_subfolders,
        z_var,
        name,
        timecol: Literal["timestep", "time", "frame"] = "timestep",
        smooth_window: int = 1,
        log_scale: bool = False,
    ):
        """Heatmap figure with one subplot per subfolder for a single variable."""
        n_rows, n_cols = self.get_subplot_layout(len(unique_subfolders))
        n_panels = n_rows * n_cols

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(6 * n_cols, 6 * n_rows),
            facecolor=self.styler.theme.figure_facecolor,
        )

        if n_panels == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes if n_cols > 1 else [axes]
        else:
            axes = axes.flatten()

        z_labels = {
            "sa_vol": "SA/Vol Ratio",
            "sa": r"Surface Area ($nm^2$)",
            "vol": r"Volume ($nm^3$)",
            "ar1": "S:M",
            "ar2": "M:L",
        }
        z_display = z_labels.get(z_var, z_var.upper())
        x_label = (
            r"$\mu$" if x_var.lower() in ("mu", "supersat", "supersaturation") else x_var.upper()
        )

        fig.suptitle(
            f"Heatmaps: {z_display} by Subfolder — {x_label} ({name.upper()})",
            fontsize=self.styler.theme.font_size_title + 2,
            fontweight=self.styler.theme.font_weight_title,
            y=0.98,
        )

        for i, subfolder in enumerate(unique_subfolders):
            if i >= len(axes):
                break
            subset = df_clean[df_clean["subfolder"] == subfolder]
            unique_x_vals = subset[x_var].nunique()
            unique_timesteps = subset[timecol].nunique()

            if unique_x_vals < 2 or unique_timesteps < 2:
                axes[i].text(
                    0.5,
                    0.5,
                    f"Insufficient data\nfor {subfolder}",
                    ha="center",
                    va="center",
                    transform=axes[i].transAxes,
                )
                axes[i].set_title(subfolder)
                continue

            try:
                noisy = subset[timecol].nunique() > subset[x_var].nunique() * smooth_window

                if noisy:
                    _c = subset[z_var]
                    if log_scale:
                        _c = np.sign(_c) * np.log10(np.abs(_c).clip(lower=1e-30))
                    hb = axes[i].hexbin(
                        subset[timecol],
                        subset[x_var],
                        C=_c,
                        gridsize=smooth_window,
                        cmap="RdYlBu_r",
                        reduce_C_function=np.mean,
                    )
                    fig.colorbar(
                        hb,
                        ax=axes[i],
                        label=f"sign\u00b7log\u2081\u2080(|{z_display}|)"
                        if log_scale
                        else z_display,
                        shrink=0.8,
                    )
                else:
                    pivot_data = subset.groupby([x_var, timecol])[z_var].mean().reset_index()
                    heatmap_data = pivot_data.pivot(index=x_var, columns=timecol, values=z_var)
                    if log_scale:
                        heatmap_data = np.sign(heatmap_data) * np.log10(
                            np.abs(heatmap_data).clip(lower=1e-30)
                        )

                    vmin = None
                    vmax = None
                    heatmap_center = heatmap_data.mean().mean()
                    if self.config.ar_limits and z_var in ["ar1", "ar2"]:
                        vmin = 0
                        vmax = 1
                        heatmap_center = None

                    log_label = (
                        f"sign\u00b7log\u2081\u2080(|{z_display}|)" if log_scale else z_display
                    )
                    sns.heatmap(
                        heatmap_data,
                        annot=False,
                        cmap="RdYlBu_r",
                        center=heatmap_center,
                        robust=True,
                        cbar_kws={"label": log_label, "shrink": 0.8},
                        linewidths=0.1,
                        linecolor="white",
                        ax=axes[i],
                        vmin=vmin,
                        vmax=vmax,
                    )

                axes[i].set_xlabel(timecol.title())
                axes[i].set_ylabel(x_label)
                axes[i].set_title(subfolder)

            except Exception:
                LOG.exception("Could not create heatmap for subfolder %s", subfolder)
                axes[i].text(
                    0.5,
                    0.5,
                    f"Error\n{subfolder}",
                    ha="center",
                    va="center",
                    transform=axes[i].transAxes,
                )
                axes[i].set_title(f"{subfolder} (Error)")

        for i in range(len(unique_subfolders), len(axes)):
            axes[i].set_visible(False)

        plt.tight_layout()

        save_name = (
            f"heatmaps_{z_var}_{x_var}_{name}_by_subfolder"
            if name
            else f"heatmaps_{z_var}_{x_var}_by_subfolder"
        )
        save_path = self.config.save_folder / f"{save_name}.png"
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

    def _create_heatmap_figure(
        self,
        df_clean,
        x_cols,
        z_var,
        name,
        mode,
        timecol: Literal["timestep", "time", "frame"] = "timestep",
        smooth_window: int = 1,
        log_scale: bool = False,
    ):
        """Helper function to create a single heatmap figure for one z variable"""

        # Set up subplot layout
        n_rows, n_cols = self.get_subplot_layout(len(x_cols))
        n_vars = n_rows * n_cols

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(6 * n_cols, 6 * n_rows),
            facecolor=self.styler.theme.figure_facecolor,
        )

        # Handle single subplot case
        if n_vars == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = [axes] if n_cols == 1 else axes
        else:
            axes = axes.flatten()

        # Set up labels
        z_labels = {
            "sa_vol": "SA/Vol Ratio",
            "sa": r"Surface Area ($nm^2$)",
            "vol": r"Volume ($nm^3$)",
            "ar1": "S:M",
            "ar2": "M:L",
        }

        # Create title
        z_display = z_labels.get(z_var, z_var.upper())
        fig.suptitle(
            f"Heatmaps: {mode.title()} Variables vs {z_display} - {name.upper()}",
            fontsize=self.styler.theme.font_size_title + 2,
            fontweight=self.styler.theme.font_weight_title,
            y=0.98,
        )

        # Create heatmaps for each x variable
        for i, x_var in enumerate(x_cols):
            if i >= len(axes):
                break

            # Check if we have enough data points
            unique_x_vals = df_clean[x_var].nunique()
            unique_timesteps = df_clean[timecol].nunique()

            if unique_x_vals < 2 or unique_timesteps < 2:
                axes[i].text(
                    0.5,
                    0.5,
                    f"Insufficient data\nfor {x_var}",
                    ha="center",
                    va="center",
                    transform=axes[i].transAxes,
                )
                axes[i].set_title(f"{x_var.upper()}")
                continue

            try:
                noisy = df_clean[timecol].nunique() > df_clean[x_var].nunique() * smooth_window

                if noisy:
                    _c = df_clean[z_var]
                    if log_scale:
                        _c = np.sign(_c) * np.log10(np.abs(_c).clip(lower=1e-30))
                    hb = axes[i].hexbin(
                        df_clean[timecol],
                        df_clean[x_var],
                        C=_c,
                        gridsize=smooth_window,
                        cmap="RdYlBu_r",
                        reduce_C_function=np.mean,
                    )
                    fig.colorbar(
                        hb,
                        ax=axes[i],
                        label=f"sign\u00b7log\u2081\u2080(|{z_display}|)"
                        if log_scale
                        else z_display,
                        shrink=0.8,
                    )
                else:
                    pivot_data = df_clean.groupby([x_var, timecol])[z_var].mean().reset_index()
                    heatmap_data = pivot_data.pivot(index=x_var, columns=timecol, values=z_var)
                    if log_scale:
                        heatmap_data = np.sign(heatmap_data) * np.log10(
                            np.abs(heatmap_data).clip(lower=1e-30)
                        )

                    vmin = None
                    vmax = None
                    heatmap_center = heatmap_data.mean().mean()
                    if self.config.ar_limits and z_var in ["ar1", "ar2"]:
                        vmin = 0
                        vmax = 1
                        heatmap_center = None

                    log_label = (
                        f"sign\u00b7log\u2081\u2080(|{z_display}|)" if log_scale else z_display
                    )
                    sns.heatmap(
                        heatmap_data,
                        annot=False,
                        cmap="RdYlBu_r",
                        center=heatmap_center,
                        robust=True,
                        cbar_kws={"label": log_label, "shrink": 0.8},
                        linewidths=0.1,
                        linecolor="white",
                        ax=axes[i],
                        vmin=vmin,
                        vmax=vmax,
                    )

                axes[i].set_xlabel(timecol.title())
                axes[i].set_ylabel(
                    r"$\mu$"
                    if x_var.lower() in ("mu", "supersat", "supersaturation")
                    else x_var.upper()
                )

            except Exception:
                LOG.exception("Could not create heatmap for %s", x_var)
                axes[i].text(
                    0.5,
                    0.5,
                    f"Error creating\nheatmap for {x_var}",
                    ha="center",
                    va="center",
                    transform=axes[i].transAxes,
                )
                axes[i].set_title(f"{x_var.upper()} (Error)")

        # Hide unused subplots
        for i in range(len(x_cols), len(axes)):
            axes[i].set_visible(False)

        plt.tight_layout()

        # Save the plot
        save_name = f"heatmaps_{z_var}_{name}" if name else f"heatmaps_{z_var}"
        save_path = self.config.save_folder / f"{save_name}.png"
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

    def plot_var_lineplots(
        self,
        df: pd.DataFrame,
        name: str = "",
        mode="parameter",
        z_vars=None,
        timecol: Literal["timestep", "time", "frame"] = "timestep",
        line_filter: dict[str, list] | None = None,
        smooth_window: int = 1,
        log_scale: bool = False,
    ):
        """
        Line plot subplots: z variable vs timestep, with each parameter value as a separate line.

        Args:
            df: DataFrame with data
            name: Name for saving
            mode: "parameter" or "energy"
            z_vars: List of y-axis variables, or None to use default ["sa_vol"]
        """
        df_clean = self._clean_columns(df)

        if z_vars is None:
            z_vars = ["sa_vol"]

        if mode == "energy":
            x_cols = [col for col in df_clean.columns if is_energy_column(col)]
            if not x_cols:
                LOG.warning("No energy columns found")
                return
        elif mode == "parameter":
            common_cols = [
                "solvent",
                "shape_name",
                "shape_file",
                "name",
                "ar1",
                "ar2",
                "sa",
                "vol",
                "sa_vol",
                "timestep",
                "time",
                "frame",
                "subfolder",
                "group",
            ]
            int_cols = [col for col in df_clean.columns if is_energy_column(col)]
            exclude = set(common_cols + int_cols)
            x_cols = [col for col in df_clean.columns if col not in exclude]
        else:
            LOG.error("Unknown mode '%s'. Use 'energy' or 'parameter'.", mode)
            return

        if not x_cols:
            LOG.warning("No x columns found for mode=%s", mode)
            return

        missing_z_vars = [z for z in z_vars if z not in df_clean.columns]
        if missing_z_vars:
            LOG.error("Missing z variables in dataframe: %s", missing_z_vars)
            return

        has_subfolders = "subfolder" in df_clean.columns
        for z_var in z_vars:
            LOG.info("Plotting Line Plot for %s Variables vs %s (%s)", mode.title(), z_var, name)
            if has_subfolders:
                unique_subfolders = sorted(df_clean["subfolder"].dropna().unique())
                for x_var in x_cols:
                    self._create_subfolder_lineplot_figure(
                        df_clean,
                        x_var,
                        unique_subfolders,
                        z_var,
                        name,
                        mode,
                        timecol,
                        line_filter,
                        smooth_window,
                        log_scale=log_scale,
                    )
            else:
                self._create_lineplot_figure(
                    df_clean,
                    x_cols,
                    z_var,
                    name,
                    mode,
                    timecol,
                    line_filter,
                    smooth_window,
                    log_scale=log_scale,
                )

    def _create_subfolder_lineplot_figure(
        self,
        df_clean,
        x_var,
        unique_subfolders,
        z_var,
        name,
        mode,
        timecol: Literal["timestep", "time", "frame"] = "timestep",
        line_filter: dict[str, list] | None = None,
        smooth_window: int = 1,
        log_scale: bool = False,
    ):
        """Single line plot: color = variable value (colorbar), marker shape = subfolder (legend)."""
        z_labels = {
            "sa_vol": "SA/Vol Ratio",
            "sa": r"Surface Area ($nm^2$)",
            "vol": r"Volume ($nm^3$)",
            "ar1": "S:M",
            "ar2": "M:L",
        }
        z_display = z_labels.get(z_var, z_var.upper())
        x_label = (
            r"$\mu$" if x_var.lower() in ("mu", "supersat", "supersaturation") else x_var.upper()
        )

        # Collect all unique variable values across all subfolders
        all_x_vals = sorted(df_clean[x_var].dropna().unique())
        if line_filter and x_var in line_filter:
            allowed = set(line_filter[x_var])
            all_x_vals = [v for v in all_x_vals if v in allowed]

        if not all_x_vals:
            LOG.warning("No values to plot for %s", x_var)
            return

        single_val = len(all_x_vals) == 1
        markers = ["o", "s", "^", "D", "v", "<", ">", "p", "*", "h"]

        if single_val:
            # One variable value: colour lines by subfolder, no colorbar
            sf_cmap = mcolors.ListedColormap(
                [
                    matplotlib.colormaps["tab10"](i / max(len(unique_subfolders) - 1, 1))
                    for i in range(len(unique_subfolders))
                ]
            )
            sf_colors = [sf_cmap(i) for i in range(len(unique_subfolders))]
            color_for_sf = dict(zip(unique_subfolders, sf_colors))

            fig, ax = plt.subplots(figsize=(9, 6), facecolor=self.styler.theme.figure_facecolor)

            for sf_idx, subfolder in enumerate(unique_subfolders):
                subset = df_clean[df_clean["subfolder"] == subfolder]
                grp = subset[subset[x_var] == all_x_vals[0]]
                if grp.empty or grp[timecol].nunique() < 2:
                    continue
                grouped = grp.groupby(timecol)[z_var].mean().reset_index()
                ax.scatter(
                    grouped[timecol],
                    grouped[z_var],
                    color=color_for_sf[subfolder],
                    marker=markers[sf_idx % len(markers)],
                    s=16,
                    alpha=0.5,
                    label=subfolder,
                )
                if smooth_window > 1:
                    y = grouped[z_var].rolling(smooth_window, center=True, min_periods=1).mean()
                    ax.plot(grouped[timecol], y, color=color_for_sf[subfolder], linewidth=1.5)

            ax.legend(fontsize=8, loc="best", framealpha=0.7)

        else:
            # Multiple variable values: colour by variable (colorbar), marker by subfolder
            all_numeric = all(isinstance(v, (int, float)) for v in all_x_vals)
            if all_numeric:
                val_arr = np.array(all_x_vals, dtype=float)
                norm = mcolors.Normalize(vmin=val_arr.min(), vmax=val_arr.max())
                cmap = matplotlib.colormaps["viridis"]

                def color_for(val):
                    return cmap(norm(float(val)))
            else:
                cmap = matplotlib.colormaps["tab10"]
                val_index = {v: i for i, v in enumerate(all_x_vals)}
                norm = mcolors.BoundaryNorm(
                    boundaries=np.arange(-0.5, len(all_x_vals) + 0.5),
                    ncolors=len(all_x_vals),
                )

                def color_for(val):
                    return cmap(val_index[val] / max(len(all_x_vals) - 1, 1))

            fig, ax = plt.subplots(figsize=(9, 6), facecolor=self.styler.theme.figure_facecolor)

            for sf_idx, subfolder in enumerate(unique_subfolders):
                subset = df_clean[df_clean["subfolder"] == subfolder]
                marker = markers[sf_idx % len(markers)]
                for val in all_x_vals:
                    grp = subset[subset[x_var] == val]
                    if grp.empty or grp[timecol].nunique() < 2:
                        continue
                    grouped = grp.groupby(timecol)[z_var].mean().reset_index()
                    ax.scatter(
                        grouped[timecol],
                        grouped[z_var],
                        color=color_for(val),
                        marker=marker,
                        s=16,
                        alpha=0.5,
                    )
                    if smooth_window > 1:
                        y = grouped[z_var].rolling(smooth_window, center=True, min_periods=1).mean()
                        ax.plot(grouped[timecol], y, color=color_for(val), linewidth=1.5)

            # Colorbar for variable values
            sm = cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.02)
            cbar.set_label(x_label)
            if not all_numeric:
                cbar.set_ticks(np.arange(len(all_x_vals)))
                cbar.set_ticklabels([str(v) for v in all_x_vals])

            # Marker legend for subfolders
            legend_handles = [
                Line2D(
                    [0],
                    [0],
                    color="grey",
                    marker=markers[i % len(markers)],
                    linestyle="-",
                    markersize=5,
                    label=sf,
                )
                for i, sf in enumerate(unique_subfolders)
            ]
            ax.legend(handles=legend_handles, fontsize=8, loc="best", framealpha=0.7)

        if log_scale:
            ax.set_yscale("symlog")

        ax.set_xlabel(timecol.title())
        ax.set_ylabel(z_display)
        ax.set_title(
            f"{z_display} vs {timecol.title()} — {x_label} ({name.upper()})",
            fontsize=self.styler.theme.font_size_title,
            fontweight=self.styler.theme.font_weight_title,
        )

        plt.tight_layout()

        save_name = (
            f"lineplots_{z_var}_{x_var}_{name}_by_subfolder"
            if name
            else f"lineplots_{z_var}_{x_var}_by_subfolder"
        )
        save_path = self.config.save_folder / f"{save_name}.png"
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

    def _create_lineplot_figure(
        self,
        df_clean,
        x_cols,
        z_var,
        name,
        mode,
        timecol: Literal["timestep", "time", "frame"] = "timestep",
        line_filter: dict[str, list] | None = None,
        smooth_window: int = 1,
        log_scale: bool = False,
    ):
        """Helper to create a line plot figure: z_var vs timecol, one line per unique x_var value."""
        n_rows, n_cols = self.get_subplot_layout(len(x_cols))
        n_vars = n_rows * n_cols

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(6 * n_cols, 5 * n_rows),
            facecolor=self.styler.theme.figure_facecolor,
        )

        if n_vars == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = [axes] if n_cols == 1 else axes
        else:
            axes = axes.flatten()

        z_labels = {
            "sa_vol": "SA/Vol Ratio",
            "sa": r"Surface Area ($nm^2$)",
            "vol": r"Volume ($nm^3$)",
            "ar1": "S:M",
            "ar2": "M:L",
        }

        z_display = z_labels.get(z_var, z_var.upper())
        fig.suptitle(
            f"Line Plots: {z_display} vs {timecol.title()} - {name.upper()}",
            fontsize=self.styler.theme.font_size_title + 2,
            fontweight=self.styler.theme.font_weight_title,
            y=0.98,
        )

        for i, x_var in enumerate(x_cols):
            if i >= len(axes):
                break

            unique_x_vals = sorted(df_clean[x_var].dropna().unique())
            if line_filter and x_var in line_filter:
                allowed = set(line_filter[x_var])
                unique_x_vals = [v for v in unique_x_vals if v in allowed]
            unique_timesteps = df_clean[timecol].nunique()

            if len(unique_x_vals) < 1 or unique_timesteps < 2:
                axes[i].text(
                    0.5,
                    0.5,
                    f"Insufficient data\nfor {x_var}",
                    ha="center",
                    va="center",
                    transform=axes[i].transAxes,
                )
                if len(x_cols) > 1:
                    axes[i].set_title(f"{x_var.upper()}")
                continue

            try:
                x_label = (
                    r"$\mu$"
                    if x_var.lower() in ("mu", "supersat", "supersaturation")
                    else x_var.upper()
                )
                use_colorbar = not (line_filter and x_var in line_filter)

                if use_colorbar:
                    all_numeric = all(isinstance(v, (int, float)) for v in unique_x_vals)
                    if all_numeric:
                        val_arr = np.array(unique_x_vals, dtype=float)
                        norm = mcolors.Normalize(vmin=val_arr.min(), vmax=val_arr.max())
                        cmap = matplotlib.colormaps["viridis"]

                        def color_for(val, cmap=cmap, norm=norm):
                            return cmap(norm(float(val)))
                    else:
                        cmap = matplotlib.colormaps["tab10"]
                        val_index = {v: idx for idx, v in enumerate(unique_x_vals)}
                        norm = mcolors.BoundaryNorm(
                            boundaries=np.arange(-0.5, len(unique_x_vals) + 0.5),
                            ncolors=len(unique_x_vals),
                        )

                        def color_for(val, cmap=cmap, val_index=val_index, x_vals=unique_x_vals):
                            return cmap(val_index[val] / max(len(x_vals) - 1, 1))
                else:
                    colors = self.styler.get_color_palette(len(unique_x_vals))

                for j, val in enumerate(unique_x_vals):
                    subset = df_clean[df_clean[x_var] == val]
                    grouped = subset.groupby(timecol)[z_var].mean().reset_index()
                    color = color_for(val) if use_colorbar else colors[j]
                    val_str = f"{val:.3g}" if isinstance(val, (int, float)) else str(val)
                    axes[i].scatter(
                        grouped[timecol],
                        grouped[z_var],
                        color=color,
                        s=9,
                        alpha=0.5,
                        **({"label": f"{x_label}={val_str}"} if not use_colorbar else {}),
                    )
                    if smooth_window > 1:
                        y = grouped[z_var].rolling(smooth_window, center=True, min_periods=1).mean()
                        axes[i].plot(grouped[timecol], y, color=color, linewidth=1.5)

                if log_scale:
                    axes[i].set_yscale("symlog")

                axes[i].set_xlabel(timecol.title())
                axes[i].set_ylabel(z_display)
                if len(x_cols) > 1:
                    axes[i].set_title(x_label)

                if use_colorbar:
                    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
                    sm.set_array([])
                    cbar = fig.colorbar(sm, ax=axes[i], pad=0.02)
                    cbar.set_label(x_label)
                    if not all_numeric:
                        cbar.set_ticks(np.arange(len(unique_x_vals)))
                        cbar.set_ticklabels([str(v) for v in unique_x_vals])
                else:
                    axes[i].legend(
                        fontsize=7,
                        loc="best",
                        framealpha=0.7,
                        title=x_label,
                    )

            except Exception:
                LOG.exception("Could not create line plot for %s", x_var)
                axes[i].text(
                    0.5,
                    0.5,
                    f"Error creating\nline plot for {x_var}",
                    ha="center",
                    va="center",
                    transform=axes[i].transAxes,
                )
                if len(x_cols) > 1:
                    axes[i].set_title(f"{x_var.upper()} (Error)")

        for i in range(len(x_cols), len(axes)):
            axes[i].set_visible(False)

        plt.tight_layout()

        save_name = f"lineplots_{z_var}_{name}" if name else f"lineplots_{z_var}"
        save_path = self.config.save_folder / f"{save_name}.png"
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

    def create_3d_zingg_plot(
        self,
        df: pd.DataFrame,
        name: str = "",
        color_var: str = "mu",
        interactive: bool = True,
        mode: str = "parameter",
        timecol: Literal["timestep", "time", "frame"] = "timestep",
    ):
        """
        Create 3D Zingg plot with timestep as z-axis and color based on specified variable

        Args:
            df: DataFrame with data
            name: Name for saving
            color_var: Variable to use for coloring ('mu', 'rand', or any other column)
            interactive: Whether to create interactive version
            mode: "parameter" or "energy" - affects how color_var is interpreted
        """
        # Clean column names
        df_clean = self._clean_columns(df)

        # Determine aspect ratio columns
        if name.startswith("cda"):
            ar_cols = [col for col in df_clean.columns if col.startswith("AspectRatio")]
            if len(ar_cols) < 2:
                LOG.error("Not enough aspect ratio columns for CDA plot")
                return
            x_col, y_col = ar_cols[0], ar_cols[1]
            x_label, y_label = ar_cols[0], ar_cols[1]
        else:
            x_col, y_col = "ar1", "ar2"
            x_label, y_label = "S:M", "M:L"

        # Check for required columns
        required_cols = [x_col, y_col, timecol]
        missing_cols = [col for col in required_cols if col not in df_clean.columns]
        if missing_cols:
            LOG.error("Missing required columns: %s", missing_cols)
            return

        # Find the color variable based on mode
        color_col = None
        if mode == "energy":
            # For energy mode, use the color_var directly (should be full column name)
            if color_var in df_clean.columns:
                color_col = color_var
        elif mode == "parameter":
            # For parameter mode, look for x_mu, x_rand, etc.
            possible_names = [f"x_{color_var}", color_var, f"{color_var}_var"]
            for possible_name in possible_names:
                if possible_name in df_clean.columns:
                    color_col = possible_name
                    break

        if color_col is None:
            LOG.warning(
                "Color variable '%s' not found in %s mode. Using timestep.",
                color_var,
                mode,
            )
            color_col = timecol  # Fallback to timestep

        # Create 3D plot
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection="3d")

        # Extract data
        x = df_clean[x_col]
        y = df_clean[y_col]
        z = df_clean[timecol]
        colors = df_clean[color_col]

        # Create scatter plot
        scatter = ax.scatter(
            x, y, z, c=colors, cmap="viridis", s=60, alpha=0.7, edgecolors="black", linewidth=0.5
        )

        # Add colorbar with appropriate label
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, aspect=30)
        if mode == "energy":
            cbar_label = "Energy (kJ/mol)" if color_is_energy_column(col) else color_col
        else:
            cbar_label = color_col.upper().replace("X_", "").replace("_", " ")
        if cbar_label.lower() in ("mu", "supersat", "supersaturation"):
            cbar_label = r"$\mu$"
        cbar.set_label(cbar_label, fontsize=12, fontweight="bold")

        # Set labels and title
        ax.set_xlabel(f"{x_label}", fontsize=12, fontweight="bold")
        ax.set_ylabel(f"{y_label}", fontsize=12, fontweight="bold")
        ax.set_zlabel(timecol.title(), fontsize=12, fontweight="bold", labelpad=0)

        title = f"3D Zingg Evolution Plot - {name.upper()}\nColored by {cbar_label} ({mode.title()} Mode)"
        ax.set_title(title, fontsize=14, fontweight="bold", pad=10)

        # Add Zingg classification regions (projected onto bottom plane)
        self._add_zingg_regions_3d(ax, x, y, z)
        self._apply_ar_limits(ax, x_is_ar=True, y_is_ar=True)

        # 3D plot appearance
        ax.grid(False)  # , alpha=0.3)
        ax.view_init(elev=20, azim=45)

        plt.tight_layout(pad=3.0)

        # Save plot
        save_name = (
            f"zingg_3d_{mode}_{color_var}_{name}" if name else f"zingg_3d_{mode}_{color_var}"
        )
        save_path = self.config.save_folder / f"{save_name}.png"
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

        # Create interactive version if requested
        if interactive:
            self._create_interactive_3d_zingg(
                df_clean, x_col, y_col, color_col, name, mode, timecol=timecol
            )

    def _add_zingg_regions_3d(self, ax, x, y, z):
        """Add Zingg classification regions to 3D plot (projected on bottom)"""
        # Get the z range for projection
        z_min = z.min()

        # Create meshgrid for the classification regions
        x_range = np.linspace(x.min(), x.max(), 50)
        y_range = np.linspace(y.min(), y.max(), 50)
        X, _ = np.meshgrid(x_range, y_range)
        _ = np.full_like(X, z_min)  # Project onto bottom plane

        # Define Zingg boundaries
        # Oblate: ar1 > 2/3, ar2 > 2/3
        # Prolate: ar1 < 2/3, ar2 < 2/3
        # Triaxial: ar1 > 2/3, ar2 < 2/3 or ar1 < 2/3, ar2 > 2/3

        # Add boundary lines projected onto bottom
        boundary_val = 2 / 3

        # Vertical line at ar1 = 2/3
        y_line = np.linspace(y.min(), y.max(), 100)
        x_line = np.full_like(y_line, boundary_val)
        z_line = np.full_like(y_line, z_min)
        ax.plot(x_line, y_line, z_line, "r--", alpha=0.5, linewidth=2, label="Zingg Boundaries")

        # Horizontal line at ar2 = 2/3
        x_line = np.linspace(x.min(), x.max(), 100)
        y_line = np.full_like(x_line, boundary_val)
        z_line = np.full_like(x_line, z_min)
        ax.plot(x_line, y_line, z_line, "r--", alpha=0.5, linewidth=2)

    def _create_interactive_3d_zingg(
        self,
        df_clean,
        x_col,
        y_col,
        color_col,
        name,
        mode="parameter",
        timecol: Literal["timestep", "time", "frame"] = "timestep",
    ):
        """Create interactive 3D Zingg plot using plotly"""
        try:
            # Create hover text
            hover_text = []
            for i in range(len(df_clean)):
                text = f"Shape: {df_clean.iloc[i].get('shape_name', 'N/A')}<br>"
                text += f"{x_col}: {df_clean.iloc[i][x_col]:.3f}<br>"
                text += f"{y_col}: {df_clean.iloc[i][y_col]:.3f}<br>"
                text += f"{timecol.title()}: {df_clean.iloc[i][timecol]}<br>"
                text += f"{color_col}: {df_clean.iloc[i][color_col]:.3f}"
                hover_text.append(text)

            # Set colorbar title based on mode
            if mode == "energy":
                colorbar_title = "Energy (kJ/mol)" if color_is_energy_column(col) else color_col
            else:
                colorbar_title = color_col.upper().replace("X_", "").replace("_", " ")

            if colorbar_title.lower() == "mu":
                colorbar_title = "μ"
            elif colorbar_title.lower() == "supersat":
                colorbar_title = "σ"

            # Create 3D scatter plot
            fig = go.Figure(
                data=go.Scatter3d(
                    x=df_clean[x_col],
                    y=df_clean[y_col],
                    z=df_clean[timecol],
                    mode="markers",
                    marker={
                        "size": 8,
                        "color": df_clean[color_col],
                        "colorscale": "Viridis",
                        "colorbar": {
                            "title": {"text": colorbar_title, "side": "right"},
                            "x": 1.05,  # push right to avoid overlap with legend
                            "y": 0.5,  # vertical center
                            "len": 0.75,  # shrink length
                        },
                        "line": {"width": 1, "color": "black"},
                    },
                    text=hover_text,
                    hovertemplate="%{text}<extra></extra>",
                    name="Data Points",
                )
            )

            # Add Zingg boundary lines
            boundary_val = 2 / 3

            # Vertical boundary line
            fig.add_trace(
                go.Scatter3d(
                    x=[boundary_val, boundary_val],
                    y=[df_clean[y_col].min(), df_clean[y_col].max()],
                    z=[df_clean[timecol].min(), df_clean[timecol].min()],
                    mode="lines",
                    line={"color": "red", "width": 4, "dash": "dash"},
                    name="Zingg Boundary (Vertical)",
                    showlegend=True,
                )
            )

            # Horizontal boundary line
            fig.add_trace(
                go.Scatter3d(
                    x=[df_clean[x_col].min(), df_clean[x_col].max()],
                    y=[boundary_val, boundary_val],
                    z=[df_clean[timecol].min(), df_clean[timecol].min()],
                    mode="lines",
                    line={"color": "red", "width": 4, "dash": "dash"},
                    name="Zingg Boundary (Horizontal)",
                    showlegend=True,
                )
            )

            # Update layout
            x_label = "S:M" if x_col == "ar1" else x_col
            y_label = "M:L" if y_col == "ar2" else y_col

            # Build scene dict
            scene_dict = {
                "xaxis_title": x_label,
                "yaxis_title": y_label,
                "zaxis_title": timecol.title(),
                "camera": {"eye": {"x": 1.5, "y": 1.5, "z": 1.2}},
            }

            if self.config.ar_limits:
                scene_dict["xaxis"] = {"title": x_label, "range": [0, 1]}
                scene_dict["yaxis"] = {"title": y_label, "range": [0, 1]}

            fig.update_layout(
                title=f"Interactive 3D Zingg Evolution - {name.upper()} ({mode.title()} Mode)",
                scene=scene_dict,
                width=1000,
                height=800,
                legend={
                    "x": 0.02,
                    "y": 0.98,  # top-left corner
                    "bgcolor": "rgba(255,255,255,0.6)",
                    "bordercolor": "black",
                    "borderwidth": 1,
                },
            )

            # Save interactive plot
            save_name = (
                f"zingg_3d_interactive_{mode}_{color_col}_{name}"
                if name
                else f"zingg_3d_interactive_{mode}_{color_col}"
            )
            save_path = self.config.save_folder / f"{save_name}.html"
            fig.write_html(save_path)

        except ImportError:
            LOG.warning("Plotly not available. Skipping interactive 3D plot.")

    def create_multiple_3d_zingg_plots(
        self,
        df: pd.DataFrame,
        name: str = "",
        mode="parameter",
        timecol: Literal["timestep", "time", "frame"] = "timestep",
    ):
        """Create 3D Zingg plots for different color variables based on mode"""

        # Clean column names first
        df_clean = self._clean_columns(df)

        # Find color variables based on mode
        if mode == "energy":
            # Energy columns
            color_cols = [col for col in df_clean.columns if is_energy_column(col)]
            if not color_cols:
                LOG.warning("No energy columns found for 3D Zingg plots")
                return
            # Use column names directly for energy mode
            color_vars = color_cols
        elif mode == "parameter":
            # Exclude common/energy columns, use the rest
            common_cols = [
                "solvent",
                "shape_name",
                "shape_file",
                "name",
                "ar1",
                "ar2",
                "sa",
                "vol",
                "sa_vol",
                "timestep",
                "frame",
            ]
            int_cols = [col for col in df_clean.columns if is_energy_column(col)]
            exclude = set(common_cols + int_cols)
            color_cols = [
                col
                for col in df_clean.columns
                if col not in exclude and pd.api.types.is_numeric_dtype(df_clean[col])
            ]

            if not color_cols:
                LOG.warning("No parameter columns found for 3D Zingg plots")
                return
            # For parameter mode, extract variable names (remove x_ prefix if present)
            color_vars = []
            for col in color_cols:
                if col.startswith("x_"):
                    var_name = col.replace("x_", "")
                else:
                    var_name = col
                color_vars.append(var_name)
        else:
            LOG.error("Unknown mode '%s'. Use 'energy' or 'parameter'.", mode)
            return

        # Create plots for each color variable
        for i, color_var in enumerate(color_vars):
            LOG.info("Creating 3D Zingg plot colored by %s (%s mode)", color_var, mode)

            # For energy mode, pass the full column name
            if mode == "energy":
                actual_color_var = color_cols[i]  # Use full column name like "Int_..."
            else:
                actual_color_var = color_var  # Use extracted name like "mu", "rand"

            self.create_3d_zingg_plot(
                df_clean,
                name=f"{name}_{mode}",
                color_var=actual_color_var,
                interactive=True,
                timecol=timecol,
            )

    def create_colored_zingg(
        self,
        df: pd.DataFrame,
        name: str = "",
        mode: Literal["energy", "parameter"] = "energy",
    ):
        """Create Zingg plots colored either by interaction energies or other parameters with consistent styling"""

        # Clean column names first
        df_clean = self._clean_columns(df)

        # Aspect ratio columns
        if name == "cda":
            ar_cols = [col for col in df_clean.columns if col.startswith("AspectRatio")]
        else:
            ar_cols = ["ar1", "ar2"]

        if len(ar_cols) < 2:
            LOG.error("Insufficient aspect ratio columns for Zingg plot")
            return

        if mode == "energy":
            # Energy columns
            color_cols = [col for col in df_clean.columns if is_energy_column(col)]
            cbar_label = "Energy (kJ/mol)"
        elif mode == "parameter":
            # Exclude common/energy columns, use the rest
            common_cols = [
                "solvent",
                "shape_name",
                "shape_file",
                "name",
                "ar1",
                "ar2",
                "sa",
                "vol",
                "sa_vol",
            ]
            int_cols = [col for col in df_clean.columns if is_energy_column(col)]
            exclude = set(common_cols + int_cols)
            color_cols = [col for col in df_clean.columns if col not in exclude]
            cbar_label = "Parameter Value"
        else:
            LOG.error("Unknown mode '%s'. Use 'energy' or 'parameter'.", mode)
            return

        if not color_cols:
            LOG.warning("No columns found for mode=%s", mode)
            return

        n_rows, n_cols = self.get_subplot_layout(len(color_cols))
        fig, axes = plt.subplots(
            n_rows, n_cols, figsize=(15, 12), facecolor=self.styler.theme.figure_facecolor
        )
        fig.suptitle(
            f"Zingg Plots by {mode.title()} - {name.upper()}",
            fontsize=self.styler.theme.font_size_title + 2,
            fontweight=self.styler.theme.font_weight_title,
            y=0.98,
        )

        if n_rows == 1:
            axes = [axes] if n_cols == 1 else axes
        else:
            axes = axes.flatten()

        LOG.info("Plotting %s Coloured Zingg Plots (%s)", mode.title(), name)

        # Choose appropriate colormap
        cmap = self.styler.get_color_palette(
            palette_type="continuous" if mode == "energy" else "diverging"
        )

        for i, col in enumerate(color_cols):
            if i >= len(axes):
                break

            ax = axes[i]
            data = pd.to_numeric(df_clean[col], errors="coerce")

            vmin = 0 if (mode == "parameter" and data.min() == -1) else None
            self._apply_ar_limits(ax, x_is_ar=True, y_is_ar=True)

            scatter = ax.scatter(
                df_clean[ar_cols[0]],
                df_clean[ar_cols[1]],
                c=data,
                cmap=cmap,
                vmin=vmin,
                alpha=self.styler.theme.alpha_scatter,
                s=self.styler.theme.marker_size,
            )

            # Apply Zingg styling
            self.styler.apply_zingg_style(ax, title=col.title(), show_legend=False)

            # Colorbar with consistent styling
            cbar = fig.colorbar(scatter, ax=ax, shrink=0.8)
            if mode == "energy":
                cbar.set_label(cbar_label, fontsize=self.styler.theme.font_size_label)
            else:
                current_label = "log[Solubility (g/L)]" if col == "solubility" else col.title()
                cbar.set_label(current_label, fontsize=self.styler.theme.font_size_label)

        # Hide unused axes
        for i in range(len(color_cols), len(axes)):
            axes[i].set_visible(False)

        plt.tight_layout()
        save_path = self.config.save_folder / f"zingg_{mode}_{name}.png"
        self.styler.save_figure(fig, save_path)

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

    def create_correlation_matrix(self, df: pd.DataFrame, name: str = ""):
        """Create correlation matrix heatmaps with consistent styling"""
        # Clean column names and rename standard columns
        df_clean = self._clean_columns(df)
        df_clean = df_clean.rename(columns={"ar1": "S:M", "ar2": "M:L"})

        numeric_df = df_clean.select_dtypes(include=[np.number]).drop(
            columns=["sa", "vol"], errors="ignore"
        )

        if numeric_df.empty:
            LOG.warning("No numeric columns found for correlation matrix")
            return

        LOG.info("Plotting Correlation Matrices (%s)", name)

        # Calculate correlations
        pearson_corr = numeric_df.corr(method="pearson")
        spearman_corr = numeric_df.corr(method="spearman")

        # Plot Pearson correlation
        fig, ax = self.styler.create_figure(figsize=(12, 10))

        # Use diverging colormap from styler
        cmap = self.styler.get_color_palette(palette_type="diverging")

        sns.heatmap(
            pearson_corr,
            annot=True,
            fmt=".2f",
            cmap=cmap,
            square=True,
            cbar=True,
            linewidths=0.5,
            ax=ax,
            cbar_kws={"shrink": 0.8},
        )

        title = (
            f"Pearson Correlation Matrix - {name.upper()}" if name else "Pearson Correlation Matrix"
        )
        self.styler.apply_correlation_style(ax, title=title)

        plt.tight_layout()
        save_path = self.config.save_folder / f"correlation_pearson_{name}.png"
        self.styler.save_figure(fig, save_path)

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

        # Plot Spearman correlation
        fig, ax = self.styler.create_figure(figsize=(12, 10))
        sns.heatmap(
            spearman_corr,
            annot=True,
            fmt=".2f",
            cmap=cmap,
            square=True,
            cbar=True,
            linewidths=0.5,
            ax=ax,
            cbar_kws={"shrink": 0.8},
        )

        title = (
            f"Spearman Correlation Matrix - {name.upper()}"
            if name
            else "Spearman Correlation Matrix"
        )
        self.styler.apply_correlation_style(ax, title=title)

        plt.tight_layout()
        save_path = self.config.save_folder / f"correlation_spearman_{name}.png"
        self.styler.save_figure(fig, save_path)

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

    def plot_solvent_effects(
        self, df: pd.DataFrame, name: str = "", exclude: list[str] | None = None
    ):
        """Plot how solvent parameters affect shape descriptors with consistent styling"""
        df_clean = self._clean_columns(df)

        if exclude:
            df_clean = df_clean[~df_clean["solvent"].isin(exclude)]

        # Calculate sphericity
        p, q = df_clean["ar1"], df_clean["ar2"]
        sphericity = (12.8 * (p**2 * q) ** (1 / 3)) / (
            1 + p * (1 + q) + 6 * (1 + p**2 * (1 + q**2)) ** 0.5
        )

        # Define parameters to plot
        common_cols = ["solvent", "ar1", "ar2", "sa", "vol", "sa_vol", "shape_name", "shape_file"]
        int_cols = [col for col in df_clean.columns if is_energy_column(col)]
        param_cols = [col for col in df_clean.columns if col not in common_cols + int_cols]

        if not param_cols:
            LOG.warning("No solvent parameters found for plotting")
            return

        n_rows, n_cols = self.get_subplot_layout(len(param_cols))
        fig, axes = plt.subplots(
            n_rows, n_cols, figsize=(15, 12), facecolor=self.styler.theme.figure_facecolor
        )
        fig.suptitle(
            f"Solvent Effects on Shape - {name.upper()}",
            fontsize=self.styler.theme.font_size_title + 2,
            fontweight=self.styler.theme.font_weight_title,
            y=0.98,
        )

        if n_rows == 1:
            axes = [axes] if n_cols == 1 else axes
        else:
            axes = axes.flatten()

        LOG.info("Plotting Solvent Parameter Effects on Shape(%s)", name)

        # Get consistent colors for different shape descriptors
        colors = self.styler.get_color_palette(3)

        for i, param in enumerate(param_cols):
            if i >= len(axes):
                break

            ax = axes[i]
            param_data = pd.to_numeric(df_clean[param], errors="coerce")

            # Create scatter plots for different shape descriptors
            ax.scatter(
                param_data,
                p,
                label="Flatness (S:M)",
                color=colors[0],
                alpha=self.styler.theme.alpha_scatter,
                s=self.styler.theme.marker_size * 0.7,
            )
            ax.scatter(
                param_data,
                q,
                label="Elongation (M:L)",
                color=colors[1],
                alpha=self.styler.theme.alpha_scatter,
                s=self.styler.theme.marker_size * 0.7,
            )
            ax.scatter(
                param_data,
                sphericity,
                label="Sphericity",
                color=colors[2],
                alpha=self.styler.theme.alpha_scatter,
                s=self.styler.theme.marker_size * 0.7,
            )

            ax.set_xlabel(
                param.replace("_", " ").title(), fontsize=self.styler.theme.font_size_label
            )
            ax.set_ylabel("Shape Ratio", fontsize=self.styler.theme.font_size_label)
            ax.set_ylim([0, 1])
            ax.set_title(
                param.replace("_", " ").title(), fontsize=self.styler.theme.font_size_title, pad=10
            )

            if i == 0:  # Add legend to first subplot only
                ax.legend(framealpha=0.9, loc="best")

        # Hide unused subplots
        for i in range(len(param_cols), len(axes)):
            axes[i].set_visible(False)

        plt.tight_layout()
        save_path = self.config.save_folder / f"solvent_effects_{name}.png"
        self.styler.save_figure(fig, save_path)

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

    def plot_interaction_energies(self, df: pd.DataFrame, name: str = ""):
        """Plot interaction energies against various parameters with consistent styling"""
        df_clean = self._clean_columns(df)
        energy_cols = [col for col in df_clean.columns if is_energy_column(col)]

        if not energy_cols:
            LOG.warning("No interaction energy columns found")
            return

        LOG.info("Plotting Solvent Parameter Effects on Energies (%s)", name)

        # Plot against solvent names
        if "solvent" in df_clean.columns:
            fig, ax = self.styler.create_figure(figsize=(12, 8))
            solvents = df_clean["solvent"].unique()

            # Get colors for different energy types
            colors = self.styler.get_color_palette(len(energy_cols))

            for i, energy_col in enumerate(energy_cols):
                for j, solvent in enumerate(solvents):
                    subset = df_clean[df_clean["solvent"] == solvent]
                    ax.scatter(
                        np.repeat(j, len(subset)),
                        subset[energy_col],
                        label=energy_col if j == 0 else "",
                        color=colors[i],
                        alpha=self.styler.theme.alpha_scatter,
                        s=self.styler.theme.marker_size * 0.8,
                    )

            ax.set_xticks(range(len(solvents)))
            ax.set_xticklabels(solvents, rotation=45, ha="right")
            ax.set_ylabel("Interaction Energy (kJ/mol)", fontsize=self.styler.theme.font_size_label)
            ax.legend(title="Interaction Type", framealpha=0.9, loc="best")
            ax.set_title(
                f"Interaction Energies by Solvent - {name.upper()}",
                fontsize=self.styler.theme.font_size_title,
                fontweight=self.styler.theme.font_weight_title,
                pad=15,
            )

            plt.tight_layout()
            save_path = self.config.save_folder / f"energies_by_solvent_{name}.png"
            self.styler.save_figure(fig, save_path)

            if self.config.show_plots:
                plt.show()
            else:
                plt.close()

    # Solvent descriptors from solvent.json use -1 to mean "not defined" (e.g. water's acidity)
    SOLVENT_PROPERTY_COLUMNS = (
        "n", "acidity", "basicity", "gamma", "dielectric", "aromatic", "halogen"
    )
    _COLOUR_BAR_LABELS = {
        "n": "Refractive index (n)",
        "dielectric": "Dielectric constant",
        "gamma": "Surface tension (γ)",
        "solubility": "log[Solubility (g/L)]",
    }

    # Shortcuts that stand for a group of columns in --colour-by / --highlight / --highlight-pic
    _COLUMN_GROUPS = {
        "energies": "energies",
        "energy": "energies",
        "ints": "energies",
        "properties": "properties",
        "props": "properties",
    }

    def _expand_column_groups(
        self, df_clean: pd.DataFrame, requested: list[str], name: str = "", warn: bool = True
    ) -> list[str]:
        """Replace group shortcuts in ``requested`` by the columns they stand for.

        ``energies`` becomes every energy column (``Int_*``, or ``Group_*`` with
        --group), in numeric order, and
        ``energies:N`` the N of them with the largest variance across the
        batch, most varied first. ``properties`` becomes every solvent property
        column present, and ``properties:N`` the N with the strongest trend
        against the aspect ratios (see :meth:`_trending_columns`). Other names
        are kept as given, and duplicates are dropped.
        """
        expanded = []
        for column in requested:
            match = re.fullmatch(r"\s*([a-z]+)\s*(?::\s*(\d+))?\s*", column.lower())
            group = self._COLUMN_GROUPS.get(match.group(1)) if match else None
            top = int(match.group(2)) if group and match.group(2) else None

            if group == "energies":
                members = natsorted(col for col in df_clean.columns if is_energy_column(col))
                if top is not None and members:
                    variances = {col: self._column_values(df_clean, col).var() for col in members}
                    # Constant energies (zero or undefined variance) can never be picked, and
                    # an Int_ column equal to one already picked (same symmetry group) is skipped
                    picked = []
                    for col in sorted(
                        (col for col in members if variances[col] > 0),
                        key=lambda col: -variances[col],
                    ):
                        values = self._column_values(df_clean, col)
                        if any(values.equals(self._column_values(df_clean, p)) for p in picked):
                            continue
                        picked.append(col)
                    members = picked[:top]
                    if warn and members:
                        LOG.info(
                            "Most varied energies (%s): %s",
                            name,
                            ", ".join(
                                f"{col} (std {variances[col] ** 0.5:.3g} kJ/mol)"
                                for col in members
                            ),
                        )
            elif group == "properties":
                members = [
                    col
                    for col in (*self.SOLVENT_PROPERTY_COLUMNS, "solubility")
                    if col in df_clean.columns
                ]
                if top is not None and members:
                    members = self._trending_columns(df_clean, members, top, name, warn)
            else:
                members = [column]
            if not members and warn:
                kind = group
                if top is not None:
                    kind = "varying energy" if group == "energies" else "trending property"
                LOG.warning("No %s columns found for '%s' (%s)", kind, column, name)
            expanded.extend(col for col in members if col not in expanded)
        return expanded

    # Fewest solvents with a value before a property's trend is worth ranking
    MIN_TREND_POINTS = 5

    def _trending_columns(
        self, df_clean: pd.DataFrame, columns: list[str], top: int, name: str, log: bool
    ) -> list[str]:
        """The ``top`` columns with the strongest trend against the aspect ratios.

        A column's trend is its largest absolute Spearman rank correlation with
        either aspect ratio, so any consistent rise or fall counts whatever the
        units. Columns with too few values, or constant ones, are left out.
        The full ranking is logged when ``log`` is set.
        """
        x, y, xlabel, ylabel = self._zingg_axes(df_clean, name)
        trends = {}
        for col in columns:
            values = self._column_values(df_clean, col)
            paired = values.notna() & x.notna() & y.notna()
            count = int(paired.sum())
            # Constant columns have no rank correlation (scipy warns on them)
            if count < self.MIN_TREND_POINTS or values[paired].nunique() < 2:
                continue
            rho = (values.corr(x, method="spearman"), values.corr(y, method="spearman"))
            score = np.nanmax(np.abs(rho)) if not np.isnan(rho).all() else np.nan
            if np.isfinite(score):
                trends[col] = (score, rho, count)

        ranked = sorted(trends, key=lambda col: -trends[col][0])
        picked = ranked[:top]
        if log and ranked:
            LOG.info("Aspect ratio trends (%s), Spearman rank correlation:", name)
            for col in ranked:
                _, (rho_x, rho_y), count = trends[col]
                LOG.info(
                    "  %s %-12s %s %+.2f  %s %+.2f  (n=%d)",
                    "*" if col in picked else " ",
                    col,
                    xlabel,
                    rho_x,
                    ylabel,
                    rho_y,
                    count,
                )
        skipped = [col for col in columns if col not in trends]
        if log and skipped:
            LOG.info(
                "  Not ranked (fewer than %d values, or constant): %s",
                self.MIN_TREND_POINTS,
                ", ".join(skipped),
            )
        return picked

    @staticmethod
    def _zingg_axes(df_clean: pd.DataFrame, name: str) -> tuple[pd.Series, pd.Series, str, str]:
        """Return the x/y aspect ratio values and their axis labels for a Zingg plot."""
        if name == "cda":
            ar_cols = [col for col in df_clean.columns if col.startswith("AspectRatio")]
            x, y = df_clean[ar_cols[0]], df_clean[ar_cols[1]]
            return x, y, ar_cols[0], ar_cols[1]
        return df_clean["ar1"], df_clean["ar2"], "S:M", "M:L"

    def _column_values(self, df_clean: pd.DataFrame, column: str) -> pd.Series:
        """Numeric values of ``column``, with the -1 solvent-property sentinel as NaN."""
        values = pd.to_numeric(df_clean[column], errors="coerce")
        if column in self.SOLVENT_PROPERTY_COLUMNS:
            values = values.mask(values == -1)
        return values

    @staticmethod
    def _available_columns(df_clean: pd.DataFrame, label_column: str) -> list[str]:
        """Columns with at least one numeric value, for 'no such column' messages."""
        return [
            col
            for col in df_clean.columns
            if col != label_column and pd.to_numeric(df_clean[col], errors="coerce").notna().any()
        ]

    def _highlight_columns(
        self,
        df_clean: pd.DataFrame,
        label_column: str,
        axes: tuple[pd.Series, pd.Series, str, str],
        requested: list[str],
        name: str,
        log: bool = True,
    ) -> dict[str, pd.Series]:
        """Map each column to find extremes in to its values, keyed by display name.

        An empty ``requested`` means the two aspect ratios on the plot's axes.
        ``ar1``/``ar2`` are displayed as ``S:M``/``M:L``. Unknown columns are
        skipped, with a warning when ``log`` is set.
        """
        x, y, xlabel, ylabel = axes
        if not requested:
            return {xlabel: x, ylabel: y}

        display_names = {"ar1": "S:M", "ar2": "M:L"} if name != "cda" else {}
        columns = {}
        for column in self._expand_column_groups(df_clean, requested, name, warn=log):
            resolved = self._resolve_column(df_clean, column)
            if resolved is None:
                if log:
                    LOG.warning(
                        "Cannot highlight %s extremes of '%s': no such column. Available: %s",
                        name,
                        column,
                        ", ".join(self._available_columns(df_clean, label_column)),
                    )
                continue
            columns[display_names.get(resolved, resolved)] = self._column_values(
                df_clean, resolved
            )
        return columns

    def _select_labelled_rows(
        self,
        df_clean: pd.DataFrame,
        label_column: str,
        axes: tuple[pd.Series, pd.Series, str, str],
        labels_to_show: list[str],
        highlight: list[str] | None,
        name: str,
        log_labels: bool = True,
        log_extremes: bool = True,
        log_heading: str = "Highlighted extremes",
    ) -> tuple[pd.Series, list[str], dict[int, list[str]]]:
        """Pick the rows to annotate and build their label text.

        Rows are picked by name from ``labels_to_show`` and, unless
        ``highlight`` is None, for holding the minimum or maximum of one of the
        ``highlight`` columns (the aspect ratios if it is empty). Extreme rows
        are tagged in their label, e.g. ``acetone (max S:M, min dielectric)``.
        ``log_labels`` and ``log_extremes`` control whether the label matches
        and the extremes (with any unknown-column warnings) are logged, so a
        selection shared by several plots is only reported once.
        Returns a boolean row mask, the label texts in mask order, and the
        extreme tags of each extreme row.
        """
        highlighted = df_clean[label_column].isin(labels_to_show)

        extremes = {}
        if highlight is not None:
            columns = self._highlight_columns(
                df_clean, label_column, axes, highlight, name, log_extremes
            )
            extremes = self._find_extremes(columns, warn=log_extremes)
            highlighted |= df_clean.index.isin(list(extremes))

        if log_labels and labels_to_show:
            LOG.info("Labels: %s", labels_to_show)
            missing = sorted(set(labels_to_show) - set(df_clean[label_column]))
            if missing:
                LOG.warning("No %s rows matched these labels: %s", label_column, missing)

        if log_extremes:
            if extremes:
                LOG.info("%s (%s):", log_heading, name)
                for row, tags in sorted(extremes.items()):
                    values = "  ".join(
                        f"{display}={series.iat[row]:.4g}" for display, series in columns.items()
                    )
                    LOG.info(
                        "  %-30s %s  [%s]",
                        df_clean.at[row, label_column],
                        values,
                        ", ".join(tags),
                    )

        label_texts = [
            f"{df_clean.at[row, label_column]} ({', '.join(extremes[row])})"
            if row in extremes
            else str(df_clean.at[row, label_column])
            for row in df_clean.index[highlighted]
        ]
        return highlighted, label_texts, extremes

    def _finish_zingg_figure(self, ax, axes, name: str, title: str):
        """Apply Zingg styling and zoom a single-panel Zingg plot onto its data."""
        x, y, xlabel, ylabel = axes
        self.styler.apply_zingg_style(ax, title=title, show_legend=False)

        # Override labels for CDA plots
        if name == "cda":
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

        # Must follow apply_zingg_style, which forces 0-1 limits
        self._zoom_to_data(ax, x, y)

    def _save_and_close(self, fig, filename: str, tight: bool = True):
        """Save a figure to the results folder, then show or close it."""
        if tight:
            fig.tight_layout()
        self.styler.save_figure(fig, self.config.save_folder / filename)

        if self.config.show_plots:
            plt.show()
        else:
            plt.close(fig)

    def create_labeled_zingg_plot(
        self,
        df: pd.DataFrame,
        labels_to_show: list[str] | None = None,
        name: str = "",
        highlight: list[str] | None = None,
    ):
        """Create Zingg plot with specific labels highlighted using consistent styling

        Points are highlighted if their name is in ``labels_to_show``, or, unless
        ``highlight`` is None, if they hold the minimum or maximum value of one
        of the ``highlight`` columns. An empty list means the two aspect ratios.
        """
        df_clean = self._clean_columns(df).reset_index(drop=True)
        axes = self._zingg_axes(df_clean, name)
        x, y, _, _ = axes

        LOG.info("Plotting Labelled Zingg Diagrams (%s)", name)

        theme = self.styler.theme
        label_column = "solvent" if "solvent" in df_clean.columns else "shape_name"
        highlight_color = theme.primary_colors["orange"]
        highlighted, label_texts, _ = self._select_labelled_rows(
            df_clean, label_column, axes, labels_to_show or [], highlight, name
        )

        fig, ax = self.styler.create_figure(figsize=(12, 10))

        # Background cloud: muted, so the highlighted points carry the figure
        ax.scatter(
            x[~highlighted],
            y[~highlighted],
            alpha=0.75,
            s=theme.marker_size * 4,
            color=theme.primary_colors["blue"],
            edgecolor="white",
            linewidth=0.4,
            zorder=2,
        )

        # Highlighted points, drawn on top and sized to be findable
        ax.scatter(
            x[highlighted],
            y[highlighted],
            alpha=0.95,
            s=theme.marker_size * 10,
            color=highlight_color,
            edgecolor=theme.primary_colors["red"],
            linewidth=1.2,
            zorder=5,
        )

        title = f"Labelled Zingg Plot - {name.upper()}" if name else "Labelled Zingg Plot"
        self._finish_zingg_figure(ax, axes, name, title)

        self._place_leader_labels(
            ax,
            list(zip(x[highlighted], y[highlighted])),
            label_texts,
            box_color=highlight_color,
            arrow_color=theme.primary_colors["red"],
        )

        self._save_and_close(fig, f"zingg_labeled_{name}.png")

    @staticmethod
    def _resolve_column(df_clean: pd.DataFrame, requested: str) -> str | None:
        """Match a user-supplied column name to a column of ``df_clean``.

        Matching ignores case and an ``x_`` prefix, and accepts ``energy1``,
        ``energy_1`` or ``int1`` as shorthand for ``Int_1``, and ``group1``,
        ``group_1`` or ``g1`` for ``Group_1``.
        """
        if requested in df_clean.columns:
            return requested

        lookup = {col.lower(): col for col in df_clean.columns}
        key = requested.strip().lower()
        candidates = [key, key.removeprefix("x_")]
        energy = re.fullmatch(r"(?:energy|int)[\s_-]*(\d+)", key)
        if energy:
            candidates.append(f"int_{energy.group(1)}")
        # Shorthands are only tried after the name itself, so e.g. "group" or "gamma" still match
        group = re.fullmatch(r"(?:group|g)[\s_-]*(\d\w*)", key)
        if group:
            candidates.append(f"group_{group.group(1)}")

        for candidate in candidates:
            if candidate in lookup:
                return lookup[candidate]
        return None

    NO_DATA_COLOUR = "#BBBBBB"

    def _colour_scale(
        self, df_clean: pd.DataFrame, column: str, name: str, label_column: str
    ) -> tuple[str, pd.Series, object, mcolors.Normalize, str] | None:
        """Resolve ``column`` and build its colour scale for a coloured Zingg plot.

        Returns the resolved column name, its values (NaN where missing), the
        colormap, the normalisation and the colour bar label, or None (with a
        warning) if the column is unknown or has no numeric values.
        """
        resolved = self._resolve_column(df_clean, column)
        if resolved is None:
            LOG.warning(
                "Cannot colour %s Zingg plot by '%s': no such column. Available: %s",
                name,
                column,
                ", ".join(self._available_columns(df_clean, label_column)),
            )
            return None

        values = self._column_values(df_clean, resolved)
        if not values.notna().any():
            LOG.warning("Column '%s' has no numeric values, skipping coloured Zingg plot", resolved)
            return None

        cmap = matplotlib.colormaps[self.styler.theme.continuous_cmap]
        norm = mcolors.Normalize(vmin=values.min(), vmax=values.max())
        if is_energy_column(resolved):
            bar_label = f"{resolved} (kJ/mol)"
        else:
            bar_label = self._COLOUR_BAR_LABELS.get(resolved, resolved.replace("_", " ").title())
        return resolved, values, cmap, norm, bar_label

    def create_zingg_coloured_by(
        self,
        df: pd.DataFrame,
        column: str,
        name: str = "",
        labels_to_show: list[str] | None = None,
        highlight: list[str] | None = None,
    ):
        """Create a single Zingg plot with points coloured by one column.

        ``column`` can be any numeric column, e.g. an interaction energy
        (``Int_1``) or a solvent property (``dielectric``). Points with no
        value are drawn in grey. Points picked by ``labels_to_show`` or
        ``highlight`` are enlarged and labelled, as in
        :meth:`create_labeled_zingg_plot`.
        """
        df_clean = self._clean_columns(df).reset_index(drop=True)
        label_column = "solvent" if "solvent" in df_clean.columns else "shape_name"

        scale = self._colour_scale(df_clean, column, name, label_column)
        if scale is None:
            return
        resolved, values, cmap, norm, bar_label = scale
        has_value = values.notna()

        LOG.info("Plotting Zingg Plot Coloured by %s (%s)", resolved, name)

        axes = self._zingg_axes(df_clean, name)
        x, y, _, _ = axes
        theme = self.styler.theme
        highlighted, label_texts, _ = self._select_labelled_rows(
            df_clean,
            label_column,
            axes,
            labels_to_show or [],
            highlight,
            name,
            log_labels=False,
            log_extremes=False,
        )

        no_data_color = self.NO_DATA_COLOUR

        fig, ax = self.styler.create_figure(figsize=(12, 10))

        def draw(rows, size, edgecolor, linewidth, zorder):
            ax.scatter(
                x[rows & has_value],
                y[rows & has_value],
                c=values[rows & has_value],
                cmap=cmap,
                norm=norm,
                s=size,
                edgecolor=edgecolor,
                linewidth=linewidth,
                zorder=zorder,
            )
            ax.scatter(
                x[rows & ~has_value],
                y[rows & ~has_value],
                color=no_data_color,
                s=size,
                edgecolor=edgecolor,
                linewidth=linewidth,
                zorder=zorder - 1,
            )

        draw(~highlighted, theme.marker_size * 4, "white", 0.4, zorder=3)
        draw(highlighted, theme.marker_size * 10, "#222222", 1.2, zorder=6)

        # Point labels sit right of the plot, so the colour bar moves underneath when there are any
        bar_kwargs = (
            {"location": "bottom", "shrink": 0.6, "pad": 0.1, "aspect": 40}
            if highlighted.any()
            else {"shrink": 0.8}
        )
        colorbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, **bar_kwargs)
        colorbar.set_label(bar_label, fontsize=theme.font_size_label)

        display = bar_label.split(" (")[0]
        title = f"Zingg Plot Coloured by {display}" + (f" - {name.upper()}" if name else "")
        safe_column = re.sub(r"[^\w.-]", "_", resolved)
        filename = f"zingg_by_{safe_column}_{name}.png"
        self._finish_zingg_figure(ax, axes, name, title)

        if not has_value.all():
            LOG.info(
                "%d %s point(s) have no %s value and are shown in grey",
                (~has_value).sum(),
                name,
                resolved,
            )
            ax.legend(
                handles=[
                    Line2D(
                        [],
                        [],
                        marker="o",
                        linestyle="",
                        color=no_data_color,
                        label=f"no {resolved} value",
                    )
                ],
                loc="best",
                framealpha=0.9,
            )

        self._place_leader_labels(
            ax,
            list(zip(x[highlighted], y[highlighted])),
            label_texts,
            box_color="#F4F4F4",
            arrow_color="#333333",
        )

        self._save_and_close(fig, filename)

    MAX_PICTURES = 12

    def _load_pictured_shapes(
        self, df_clean: pd.DataFrame, rows: list[int], label_column: str
    ) -> dict[int, np.ndarray]:
        """Load the shape coordinates behind ``rows``, skipping any that cannot be read.

        Shapes are cached, since the plain and coloured picture plots draw the same ones.
        """
        if not hasattr(self, "_shape_cache"):
            self._shape_cache: dict[tuple[Path, int], np.ndarray] = {}
        shapes = {}
        for row in rows:
            label = df_clean.at[row, label_column]
            entry = df_clean.at[row, "shape_file"]
            if pd.isna(entry):
                LOG.warning("No shape file recorded for %s", label)
                continue
            path = self.config.resolve_shape_file(entry)
            if not path.exists():
                LOG.warning("Shape file for %s not found: %s", label, path)
                continue
            frame = int(df_clean.at[row, "frame"]) if "frame" in df_clean.columns else -1
            if (path, frame) in self._shape_cache:
                shapes[row] = self._shape_cache[(path, frame)]
                continue
            try:
                shapes[row] = self._shape_cache[(path, frame)] = load_shape_points(path, frame)
            except Exception:
                LOG.exception("Could not load the shape for %s from %s", label, path)
        return shapes

    def create_zingg_pictures(
        self,
        df: pd.DataFrame,
        labels_to_show: list[str] | None = None,
        name: str = "",
        highlight: list[str] | None = None,
        log_labels: bool = False,
        log_extremes: bool = False,
        colour_by: str | None = None,
    ):
        """Create a Zingg plot with a picture of the crystal behind each labelled point.

        Points are picked as in :meth:`create_labeled_zingg_plot`. Each one is
        joined by a leader line to a drawing of its crystal (the convex hull of
        its shape file, see :mod:`cgpytools.plot.crystal_render`) in a column
        right of the plot. Needs the ``shape_file`` column the analyses record.
        The log flags are passed to :meth:`_select_labelled_rows`; they are off
        by default since the labelled plot usually reports the same selection.
        With ``colour_by``, points are coloured by that column as in
        :meth:`create_zingg_coloured_by` and the plot is saved as
        ``zingg_pictures_by_<column>_<name>.png``.
        """
        df_clean = self._clean_columns(df).reset_index(drop=True)
        if "shape_file" not in df_clean.columns:
            LOG.warning(
                "Cannot draw crystal pictures for %s: the results have no shape_file column. %s",
                name,
                "CDA results have no crystal geometry."
                if name == "cda"
                else "Re-run the analysis from --input_dir to record the shape files.",
            )
            return

        label_column = "solvent" if "solvent" in df_clean.columns else "shape_name"
        scale = None
        if colour_by is not None:
            scale = self._colour_scale(df_clean, colour_by, name, label_column)
            if scale is None:
                return
        axes = self._zingg_axes(df_clean, name)
        x, y, _, _ = axes
        highlighted, _, extremes = self._select_labelled_rows(
            df_clean,
            label_column,
            axes,
            labels_to_show or [],
            highlight,
            name,
            log_labels=log_labels,
            log_extremes=log_extremes,
            log_heading="Pictured extremes",
        )

        rows = list(df_clean.index[highlighted])
        if not rows:
            LOG.warning("No points selected for crystal pictures (%s)", name)
            return
        if len(rows) > self.MAX_PICTURES:
            LOG.warning(
                "%d points selected for crystal pictures (%s); drawing the first %d",
                len(rows),
                name,
                self.MAX_PICTURES,
            )
            rows = rows[: self.MAX_PICTURES]

        if scale is None:
            LOG.info("Plotting Zingg Plot with Crystal Pictures (%s)", name)
        else:
            LOG.info(
                "Plotting Zingg Plot with Crystal Pictures Coloured by %s (%s)", scale[0], name
            )
        shapes = self._load_pictured_shapes(df_clean, rows, label_column)
        if not shapes:
            LOG.warning("No shape files could be loaded, skipping crystal pictures (%s)", name)
            return

        # Pictures sit in one column right of the plot, where there are no axis labels
        # for the leader lines to cross, ordered top to bottom like their points
        ordered = sorted(shapes, key=lambda row: -y.iat[row])

        theme = self.styler.theme
        # Grow the figure for many pictures; the box aspect stops the plot stretching with it
        height = max(9.0, 1.5 * len(ordered) + 1.5)
        fig = plt.figure(figsize=(12.5, height), facecolor=theme.figure_facecolor)
        grid = fig.add_gridspec(1, 2, width_ratios=[4, 1.5], wspace=0.15)
        ax = fig.add_subplot(grid[0, 0])
        ax.set_box_aspect(0.85)

        pictured = df_clean.index.isin(list(shapes))
        if scale is None:
            ax.scatter(
                x[~pictured],
                y[~pictured],
                alpha=0.75,
                s=theme.marker_size * 4,
                color=theme.primary_colors["blue"],
                edgecolor="white",
                linewidth=0.4,
                zorder=2,
            )
            ax.scatter(
                x[pictured],
                y[pictured],
                alpha=0.95,
                s=theme.marker_size * 10,
                color=theme.primary_colors["orange"],
                edgecolor=theme.primary_colors["red"],
                linewidth=1.2,
                zorder=5,
            )
            title = "Crystal Shapes on the Zingg Plot"
            filename = f"zingg_pictures_{name}.png"
            leader_colour = theme.primary_colors["red"]
        else:
            resolved, values, cmap, norm, bar_label = scale
            has_value = values.notna()
            for rows_mask, size, edge, width, zorder in (
                (~pictured, theme.marker_size * 4, "white", 0.4, 3),
                (pictured, theme.marker_size * 10, "#222222", 1.2, 6),
            ):
                ax.scatter(
                    x[rows_mask & has_value],
                    y[rows_mask & has_value],
                    c=values[rows_mask & has_value],
                    cmap=cmap,
                    norm=norm,
                    s=size,
                    edgecolor=edge,
                    linewidth=width,
                    zorder=zorder,
                )
                ax.scatter(
                    x[rows_mask & ~has_value],
                    y[rows_mask & ~has_value],
                    color=self.NO_DATA_COLOUR,
                    s=size,
                    edgecolor=edge,
                    linewidth=width,
                    zorder=zorder - 1,
                )
            colorbar = fig.colorbar(
                cm.ScalarMappable(norm=norm, cmap=cmap),
                ax=ax,
                location="bottom",
                shrink=0.6,
                pad=0.08,
                aspect=40,
            )
            colorbar.set_label(bar_label, fontsize=theme.font_size_label)
            title = f"Crystal Shapes Coloured by {bar_label.split(' (')[0]}"
            safe_column = re.sub(r"[^\w.-]", "_", resolved)
            filename = f"zingg_pictures_by_{safe_column}_{name}.png"
            leader_colour = "#333333"
        title += f" - {name.upper()}" if name else ""
        self._finish_zingg_figure(ax, axes, name, title)

        cells = grid[0, 1].subgridspec(len(ordered), 1, hspace=0.5)
        for i, row in enumerate(ordered):
            pic_ax = fig.add_subplot(cells[i, 0])
            draw_crystal(pic_ax, shapes[row], color=theme.primary_colors["blue"])

            caption = str(df_clean.at[row, label_column])
            if "frame" in df_clean.columns:
                caption += f" (frame {df_clean.at[row, 'frame']})"
            elif "timestep" in df_clean.columns:
                caption += f" (timestep {df_clean.at[row, 'timestep']})"
            if row in extremes:
                caption += "\n" + ", ".join(extremes[row])
            pic_ax.set_title(caption, fontsize=theme.font_size_tick, pad=3)

            point = (x.iat[row], y.iat[row])
            if np.isfinite(point).all():
                fig.add_artist(
                    ConnectionPatch(
                        xyA=point,
                        coordsA=ax.transData,
                        xyB=(0.0, 0.5),
                        coordsB=pic_ax.transAxes,
                        color=leader_colour,
                        linewidth=1.0,
                        alpha=0.8,
                        shrinkA=5,
                        shrinkB=4,
                        zorder=6,
                    )
                )

        # tight_layout does not handle nested grids; the margins are set directly instead
        fig.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.07)
        self._save_and_close(fig, filename, tight=False)

    def create_annotated_zingg_plots(
        self, df: pd.DataFrame, labels_to_show: list[str] | None = None, name: str = ""
    ):
        """Create the labelled, colour-by and picture Zingg plots requested in the config."""
        highlight = self.config.highlight
        labelled_plot = bool(labels_to_show) or highlight is not None
        if labelled_plot:
            self.create_labeled_zingg_plot(df, labels_to_show, name=name, highlight=highlight)
        colour_by = self._expand_column_groups(
            self._clean_columns(df), self.config.colour_by, name
        )
        for column in colour_by:
            self.create_zingg_coloured_by(
                df, column, name=name, labels_to_show=labels_to_show, highlight=highlight
            )

        pictures = self.config.highlight_pictures
        if pictures is not None:
            # Columns given to --highlight-pic win; otherwise picture what --highlight or
            # --labels picked; with neither, picture the aspect ratio extremes
            if pictures:
                picture_highlight = pictures
            elif highlight is not None or labels_to_show:
                picture_highlight = highlight
            else:
                picture_highlight = []
            self.create_zingg_pictures(
                df,
                labels_to_show,
                name=name,
                highlight=picture_highlight,
                # Report only what the labelled plot has not already reported
                log_labels=not labelled_plot,
                log_extremes=picture_highlight != highlight,
            )
            # With --colour-by too, the same pictures again with points coloured by each column
            for column in colour_by:
                self.create_zingg_pictures(
                    df, labels_to_show, name=name, highlight=picture_highlight, colour_by=column
                )

    def set_plot_style(self, style: str, custom_theme: PlotTheme = None):
        """Change the global plot style"""
        self.styler.reset_style()
        self.styler = GlobalPlotStyler(theme=custom_theme, style=style)
        LOG.info("Plot style changed to: %s", style)

    def __del__(self):
        """Cleanup: reset matplotlib settings when object is destroyed"""
        if hasattr(self, "styler"):
            self.styler.reset_style()


class CrystalShapeAnalysisPipeline:
    """Main pipeline class that orchestrates the entire analysis"""

    def __init__(
        self,
        config: ShapeAnalysisConfig,
        plot_style: Literal["modern", "classic", "minimal", "dark", "publication"] = "modern",
        custom_theme: PlotTheme = None,
    ):
        self.config = config
        self.crystal_analyser = CrystalShapeAnalyser(config)
        self.plots = ShapePlots(config, plot_style=plot_style, custom_theme=custom_theme)

    def group_energies(self, df: pd.DataFrame) -> pd.DataFrame:
        """With --group, give the plots one energy per symmetry group of interactions.

        The results CSV keeps the per-line ``Int_*`` columns; only the plotted
        frame is grouped, so results can be re-plotted either way. The dimer
        behind each group is also exported to ``key_dimers``, see
        :meth:`export_key_dimers`.
        """
        if not self.config.group_energies or df.empty:
            return df

        netfile = EnergyDataLoader.find_net_file(df, self.config)
        if netfile is None:
            LOG.warning(
                "--group: no net.txt found beside the shape files or under %s; "
                "plotting ungrouped energies",
                self.config.input_dir,
            )
            return df
        groups = EnergyDataLoader.read_energy_groups(netfile)
        if groups is None:
            return df

        if any(re.fullmatch(r"Int_\d+", col) for col in df.columns):
            df = EnergyDataLoader.group_energy_columns(df, groups, netfile)
        else:
            LOG.warning("--group: no interaction energies in these results (use --energies)")
        self.export_key_dimers(groups, df)
        return df

    def export_key_dimers(self, groups: dict[str, list[dict]], df: pd.DataFrame):
        """Copy, describe and draw the dimer behind each energy group (see :mod:`dimers`)."""
        folder = find_dimers_folder(self.config.input_dir)
        if folder is None:
            LOG.info(
                "--group: no folder ending in 'dimers' in %s, so no key dimers exported",
                self.config.input_dir,
            )
            return
        LOG.info("Exporting key dimers from %s", folder)
        energy_cols = [name for name in groups if name in df.columns]
        pairs = {name: [line["pair"] for line in lines] for name, lines in groups.items()}
        export_key_dimers(
            {name: [line["dimer"] for line in lines] for name, lines in groups.items()},
            folder,
            self.config.save_folder / "key_dimers",
            energies=df[energy_cols] if energy_cols else None,
            multiplicity={name: len(lines) for name, lines in groups.items()},
            distances={name: lines[0]["r"] for name, lines in groups.items()},
            pairs=pairs,
            crystal=self.load_crystal(folder, pairs),
        )

    def load_crystal(self, dimers_folder: Path, pairs: dict[str, list]) -> CrystalMotifs | None:
        """The crystal the key dimers come from, to give their directions in it."""
        cif = find_cif(dimers_folder, self.config.input_dir)
        if cif is None:
            LOG.info(
                "--group: no CIF beside %s, so no interaction directions for the key dimers",
                dimers_folder,
            )
            return None
        try:
            crystal = CrystalMotifs(cif, [p for lines in pairs.values() for p in lines if p])
        except Exception as error:  # noqa: BLE001 - optional extra, never stops the export
            LOG.warning(
                "--group: could not place the net file's dimers in %s (%s), so no interaction "
                "directions for the key dimers",
                cif,
                error,
            )
            return None
        LOG.info("Placing key dimers in the crystal of %s", cif)
        return crystal

    def run_general_analysis(
        self,
        shape_files: list[Path],
        energy_csv: Path | None = None,
        labels_to_show: list[str] | None = None,
    ) -> pd.DataFrame:
        """Run analysis for general shapes (not from solvent screening)"""
        LOG.info("Starting general shape analysis")

        # Analyse shapes
        df = self.crystal_analyser.analyse_general_shapes(
            shape_files,
            energy_csv,
        )

        if df.empty:
            LOG.error("No shapes were successfully analysed")
            return df

        # Create visualisations
        self.plots.create_zingg_plot(df, name="general")
        self.plots.create_correlation_matrix(df, name="general")

        # Labelled and colour-by Zingg plots, if requested
        self.plots.create_annotated_zingg_plots(df, labels_to_show, name="general")

        # If energy data is available, create energy-colored plots
        param_cols = [col for col in df.columns if col.startswith((*ENERGY_PREFIXES, "x_", "X_"))]
        if param_cols:
            self.plots.create_colored_zingg(df, name="general", mode="parameter")

        LOG.info("General shape analysis completed")
        return df

    def run_solvent_analysis(
        self,
        shape_files: list[Path],
        solvent_json: Path,
        occ_outputs: list[Path] | None = None,
        get_energies: bool = False,
        exclude_solvents: list[str] | None = None,
        labels_to_show: list[str] | None = None,
    ) -> pd.DataFrame:
        """Run analysis for solvent screening results"""
        LOG.info("Starting solvent screening analysis")

        # Analyse shapes with solvent data
        df = self.crystal_analyser.analyse_solvent_shapes(
            shape_files, solvent_json, occ_outputs, get_energies
        )

        if df.empty:
            LOG.error("No solvent shapes were successfully analysed")
            return df
        df = self.group_energies(df)

        # Create visualisations
        self.plots.create_zingg_plot(df, name="solvent")
        self.plots.create_colored_zingg(df, name="solvent", mode="parameter")
        self.plots.create_correlation_matrix(df, name="solvent")
        self.plots.plot_solvent_effects(df, name="solvent", exclude=exclude_solvents)

        # Labelled and colour-by Zingg plots, if requested
        self.plots.create_annotated_zingg_plots(df, labels_to_show, name="solvent")

        # Create interaction energy plots if available
        energy_cols = [col for col in df.columns if is_energy_column(col)]
        if energy_cols:
            self.plots.plot_interaction_energies(df, name="solvent")
            self.plots.create_colored_zingg(df, name="solvent", mode="energy")

        LOG.info("Solvent screening analysis completed")
        return df

    def run_cda_analysis(
        self,
        cda_files: list[Path],
        directions: list[str],
        get_energies: bool = False,
        labels_to_show: list[str] | None = None,
    ) -> pd.DataFrame:
        """Run analysis for CDA simulation results"""
        LOG.info("Starting CDA analysis")

        if len(directions) != 3:
            raise ValueError(
                f"Exactly 3 directions required for CDA analysis, got {len(directions)}"
            )

        # Analyse CDA results
        df = self.crystal_analyser.analyse_cda_shapes(cda_files, directions, get_energies)

        if df.empty:
            LOG.error("No CDA shapes were successfully analysed")
            return df
        df = self.group_energies(df)

        # Create visualisations
        self.plots.create_zingg_plot(df, name="cda")

        # Labelled and colour-by Zingg plots, if requested
        self.plots.create_annotated_zingg_plots(df, labels_to_show, name="cda")

        # Create energy plots if available
        energy_cols = [col for col in df.columns if is_energy_column(col)]
        if energy_cols:
            self.plots.create_colored_zingg(df, name="cda", mode="energy")

        LOG.info("CDA analysis completed")
        return df

    def run_size_wulff_analysis(
        self,
        generation_results: dict[str, dict],
        energy_csvs: list[Path] | None = None,
        relative: bool = True,
        labels_to_show: list[str] | None = None,
        line_filter: dict[str, list] | None = None,
        smooth_window: int = 1,
        timecol: Literal["timestep", "time", "frame"] = "timestep",
        log_scale: bool = False,
    ) -> pd.DataFrame:
        """
        Run analysis for size/Wulff shapes from generation results

        Args:
            generation_results: Results from generation step containing timestep data
            energy_csvs: Optional list of CSV paths with energy data (one per subfolder group)
            relative: Whether surface area and volume are zeroed at initial frame
            labels_to_show: Optional list of labels for labeled plots

        Returns:
            DataFrame with shape analysis results including growth metrics
        """
        LOG.info("Starting size/Wulff shape analysis")

        # Analyze shapes using the existing method
        df = self.crystal_analyser.analyse_size_wulff_shapes(
            generation_results,
            energy_csvs,
            relative,
        )
        LOG.info(df)

        if df.empty:
            LOG.error("No size/Wulff shapes were successfully analyzed")
            return df

        # Create basic visualisations
        self.plots.create_zingg_plot(df, name="wulff")
        self.plots.create_colored_zingg(df, name="wulff", mode="parameter")
        self.plots.create_correlation_matrix(df, name="wulff")

        # Labelled and colour-by Zingg plots, if requested
        self.plots.create_annotated_zingg_plots(df, labels_to_show, name="wulff")

        self.plots.plot_var_heatmaps(
            df,
            name="wulff",
            mode="parameter",
            z_vars=["ar1", "ar2"],
            timecol=timecol,
            smooth_window=smooth_window,
        )
        self.plots.plot_var_heatmaps(
            df,
            name="wulff",
            mode="parameter",
            z_vars=["sa", "vol"],
            timecol=timecol,
            smooth_window=smooth_window,
            log_scale=log_scale,
        )
        self.plots.plot_var_heatmaps(
            df,
            name="wulff",
            mode="parameter",
            z_vars=["sa_vol"],
            timecol=timecol,
            smooth_window=smooth_window,
        )
        self.plots.plot_var_lineplots(
            df,
            name="wulff",
            mode="parameter",
            z_vars=["ar1", "ar2"],
            timecol=timecol,
            line_filter=line_filter,
            smooth_window=smooth_window,
        )
        self.plots.plot_var_lineplots(
            df,
            name="wulff",
            mode="parameter",
            z_vars=["sa", "vol"],
            timecol=timecol,
            line_filter=line_filter,
            smooth_window=smooth_window,
            log_scale=log_scale,
        )
        self.plots.plot_var_lineplots(
            df,
            name="wulff",
            mode="parameter",
            z_vars=["sa_vol"],
            timecol=timecol,
            line_filter=line_filter,
            smooth_window=smooth_window,
        )

        self.plots.create_multiple_3d_zingg_plots(df, name="wulff", mode="parameter")

        # Check if energy data is available

        energy_cols = [col for col in df.columns if is_energy_column(col)]
        if energy_cols:
            # You could also do energy columns vs timestep with different z variables
            self.plots.create_colored_zingg(df, name="wulff", mode="energy")
            self.plots.plot_var_heatmaps(
                df,
                name="wulff",
                mode="energy",
                z_vars=["sa_vol"],
                timecol=timecol,
                smooth_window=smooth_window,
            )
            self.plots.plot_var_lineplots(
                df,
                name="wulff",
                mode="energy",
                z_vars=["sa_vol"],
                timecol=timecol,
                line_filter=line_filter,
                smooth_window=smooth_window,
            )
            self.plots.create_multiple_3d_zingg_plots(df, name="wulff", mode="energy")
        else:
            LOG.warning("No energy data available, skipping potential plots")

        LOG.info("Size/Wulff shape analysis completed for %s data points", len(df))
        return df

    def run_movie_analysis(
        self,
        movie_files: list[Path],
        energy_csv: Path | None = None,
        labels_to_show: list[str] | None = None,
        line_filter: dict[str, list] | None = None,
        smooth_window: int = 1,
    ) -> pd.DataFrame:
        """
        Run analysis for shape movies, tracking shape metrics over frames.
        """
        LOG.info("Starting movie shape analysis")
        df = self.crystal_analyser.analyse_movies(movie_files, energy_csv)

        if df.empty:
            LOG.error("No movies were successfully analysed")
            return df

        # Create basic visualisations
        self.plots.create_zingg_plot(df, name="movies")
        self.plots.create_colored_zingg(df, name="movies", mode="parameter")
        self.plots.create_correlation_matrix(df, name="movies")

        # Labelled and colour-by Zingg plots, if requested
        self.plots.create_annotated_zingg_plots(df, labels_to_show, name="movies")

        self.plots.plot_var_heatmaps(
            df, name="movies", mode="parameter", z_vars=["ar1", "ar2"], timecol="frame"
        )
        self.plots.plot_var_heatmaps(
            df, name="movies", mode="parameter", z_vars=["sa", "vol"], timecol="frame"
        )
        self.plots.plot_var_heatmaps(
            df, name="movies", mode="parameter", z_vars=["sa_vol"], timecol="frame"
        )
        self.plots.plot_var_lineplots(
            df,
            name="movies",
            mode="parameter",
            z_vars=["ar1", "ar2"],
            timecol="frame",
            line_filter=line_filter,
            smooth_window=smooth_window,
        )
        self.plots.plot_var_lineplots(
            df,
            name="movies",
            mode="parameter",
            z_vars=["sa", "vol"],
            timecol="frame",
            line_filter=line_filter,
            smooth_window=smooth_window,
        )
        self.plots.plot_var_lineplots(
            df,
            name="movies",
            mode="parameter",
            z_vars=["sa_vol"],
            timecol="frame",
            line_filter=line_filter,
            smooth_window=smooth_window,
        )

        self.plots.create_multiple_3d_zingg_plots(
            df, name="movies", mode="parameter", timecol="frame"
        )

        # Check if energy data is available

        energy_cols = [col for col in df.columns if is_energy_column(col)]
        if energy_cols:
            # You could also do energy columns vs timestep with different z variables
            self.plots.create_colored_zingg(df, name="movies", mode="energy")
            self.plots.plot_var_heatmaps(
                df, name="movies", mode="energy", z_vars=["sa_vol"], timecol="frame"
            )
            self.plots.plot_var_lineplots(
                df,
                name="movies",
                mode="energy",
                z_vars=["sa_vol"],
                timecol="frame",
                line_filter=line_filter,
                smooth_window=smooth_window,
            )
            self.plots.create_multiple_3d_zingg_plots(
                df, name="movies", mode="energy", timecol="frame"
            )
        else:
            LOG.warning("No energy data available, skipping potential plots")

        LOG.info(
            "Movie shape analysis completed for %d frames across %d movies",
            len(df),
            df["shape_name"].nunique(),
        )
        return df

    def run_wulff_analysis(
        self,
        wulff_files: list[Path],
        energy_csv: Path | None = None,
        labels_to_show: list[str] | None = None,
        line_filter: dict[str, list] | None = None,
        smooth_window: int = 1,
    ) -> pd.DataFrame:
        """Run analysis for PLY Wulff shape files (same pipeline as general XYZ analysis)"""
        LOG.info("Starting Wulff PLY shape analysis")

        # Analyse shapes using the general shapes method (works for any mesh format)
        df = self.crystal_analyser.analyse_general_shapes(
            wulff_files,
            energy_csv,
        )

        if df.empty:
            LOG.error("No Wulff shapes were successfully analysed")
            return df

        # Create visualisations
        self.plots.create_zingg_plot(df, name="wulff_ply")
        self.plots.create_correlation_matrix(df, name="wulff_ply")

        # Labelled and colour-by Zingg plots, if requested
        self.plots.create_annotated_zingg_plots(df, labels_to_show, name="wulff_ply")

        # If energy data is available, create energy-colored plots
        param_cols = [col for col in df.columns if col.startswith((*ENERGY_PREFIXES, "x_", "X_"))]
        if param_cols:
            self.plots.create_colored_zingg(df, name="wulff_ply", mode="parameter")

        LOG.info("Wulff PLY shape analysis completed")
        return df


def build_parser() -> argparse.ArgumentParser:
    """Build the command line interface for the screening tool."""
    parser = argparse.ArgumentParser(description="Crystal Shape Analysis Tool")

    # Input/Output arguments
    parser.add_argument("-i", "--input_dir", type=Path, help="Path to input folder")
    parser.add_argument("-o", "--output_dir", type=Path, help="Path to output folder")
    parser.add_argument("-r", "--results_dir", type=Path, help="Path to existing results CSV files")

    # Data source arguments
    parser.add_argument(
        "-s",
        "--solvent_json",
        type=Path,
        default="solvent.json",
        help="Path to solvent properties JSON file",
    )
    parser.add_argument(
        "-e",
        "--energy_csv",
        nargs="*",
        type=Path,
        help="Path(s) to CSV file(s) containing energy data (one per subfolder group)",
    )
    parser.add_argument(
        "-c",
        "--crystallography",
        type=Path,
        help="Path to JSON file containing crystallographic information",
    )
    parser.add_argument(
        "--energy_name_col",
        default="name",
        help="Column name in energy CSV that contains shape names",
    )

    # Analysis type flags
    parser.add_argument(
        "--general", action="store_true", help="Run general shape analysis (not solvent-specific)"
    )
    parser.add_argument("--solvent", action="store_true", help="Run solvent screening analysis")
    parser.add_argument("--wulff", action="store_true", help="Analyse PLY Wulff shape files")
    parser.add_argument("--cda", action="store_true", help="Analyse CDA results")
    parser.add_argument("--occ", action="store_true", help="Include OCC solubility data")

    # Processing options
    parser.add_argument(
        "--energies", action="store_true", help="Extract interaction energies from net files"
    )
    parser.add_argument(
        "--group",
        action="store_true",
        help="Plot one energy per symmetry group of interactions in the net files (Group_1, "
        "Group_2, ... by net file ID) instead of one per interaction line (Int_1, Int_2, ...). "
        "For --solvent and --cda; the results CSV keeps the Int_ columns",
    )
    parser.add_argument(
        "--movies", action="store_true", help="Extract growth information from XYZ movie files"
    )
    parser.add_argument(
        "--size", action="store_true", help="Extract growth information from size.csv files"
    )
    parser.add_argument(
        "--no-relative",
        action="store_true",
        help="Store absolute surface area and volume instead of values relative to the initial frame",
    )
    parser.add_argument(
        "--log-scale",
        action="store_true",
        help="Plot surface area and volume on a log scale",
    )
    parser.add_argument(
        "--wulff-interval",
        type=int,
        default=10,
        help="Generate a Wulff shape mesh every n step of the size.csv file",
    )
    parser.add_argument(
        "--box", action="store_true", help="Use bounding box to calculate Zingg plot."
    )
    parser.add_argument(
        "--lmax", type=int, default=20, help="Resolution for Wulff construction rendering"
    )

    # CDA-specific arguments
    parser.add_argument(
        "-d",
        "--directions",
        nargs=3,
        help="Three crystallographic directions for CDA analysis (S:M:L order)",
    )

    # visualisation options
    parser.add_argument("--show", action="store_true", help="Display plots interactively")
    parser.add_argument(
        "--ar-limits",
        action="store_true",
        help="Set aspect ratio axes to always show 0-1 range in plots",
    )
    parser.add_argument(
        "--zoom",
        action="store_true",
        help="Zoom the labelled, --colour-by and --highlight-pic Zingg plots onto the data "
        "instead of the full 0-1 Zingg diagram",
    )
    parser.add_argument(
        "-l", "--labels", nargs="+", default=[], help="Specific shapes/solvents to label in plots"
    )
    parser.add_argument(
        "--highlight",
        nargs="*",
        metavar="COLUMN",
        help="Label the shapes/solvents with the minimum and maximum of each named column "
        "(the aspect ratios if none are named) on the Zingg plots, and print them to the log, "
        "e.g. --highlight dielectric acidity ar1 ar2. Takes the --colour-by shortcuts "
        "(energies, energies:N, properties, properties:N). Combines with --labels",
    )
    parser.add_argument(
        "--highlight-pic",
        nargs="*",
        metavar="COLUMN",
        help="Draw a picture of the crystal behind the minimum and maximum of each named column "
        "beside the Zingg plot, e.g. --highlight-pic dielectric energies:2. With no columns, "
        "picture what --highlight/--labels pick (the aspect ratio extremes if neither is given). "
        "Points from --labels are always pictured",
    )
    parser.add_argument(
        "--colour-by",
        "--color-by",
        dest="colour_by",
        nargs="+",
        default=[],
        metavar="COLUMN",
        help="Make a Zingg plot coloured by each named column, e.g. --colour-by dielectric "
        "acidity Int_1 (energy1 also works for Int_1). 'energies' makes one per interaction "
        "energy, 'energies:N' one for each of the N energies that vary most across the batch, "
        "'properties' one per solvent property, and 'properties:N' one for each of the N "
        "properties that trend most with the aspect ratios. Points from --labels/--highlight "
        "are labelled on these plots too",
    )
    parser.add_argument(
        "--overview",
        type=int,
        metavar="N",
        help="Quick insight into a system: shorthand for --group --highlight energies:N "
        "--highlight-pic energies:N --colour-by properties:N, i.e. the N energies that vary "
        "most and the N solvent properties that trend most with shape. The crystal pictures "
        "are drawn plain and coloured by each property. Adds to any of those flags given "
        "explicitly. Use with --energies",
    )
    parser.add_argument(
        "--exclude", nargs="+", default=[], help="Solvents to exclude from analysis"
    )
    parser.add_argument(
        "--line-filter",
        nargs="+",
        default=[],
        metavar="COL=VAL1,VAL2",
        help="Filter lines in line plots by column value, e.g. --line-filter mu=0.1,0.5 ss=1.2",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=1,
        metavar="N",
        help="Averaging window size for smoothing line plots (default: 1 = no smoothing)",
    )
    parser.add_argument(
        "--font-scale",
        type=float,
        default=1.0,
        metavar="FACTOR",
        help="Scale every font size in the plots, e.g. 1.5 for 50%% larger text",
    )
    parser.add_argument(
        "--border-width",
        type=float,
        metavar="WIDTH",
        help="Width of the axes frame in points (default 1.2). "
        "The Zingg 2/3 classification lines scale with it",
    )
    parser.add_argument(
        "--border-color",
        metavar="COLOR",
        help="Colour of the axes frame (default #CCCCCC). "
        "Worth darkening, e.g. #333333, when increasing --border-width",
    )
    parser.add_argument(
        "--plot-style",
        choices=["modern", "classic", "minimal", "dark", "publication"],
        default="modern",
        help="Base plot style (default: modern)",
    )
    parser.add_argument(
        "--timecol",
        choices=["timestep", "time", "frame"],
        default="timestep",
        help="Column to use as the time axis in line/heat plots (default: timestep). "
        "Use 'time' for the real time value from size.csv.",
    )

    return parser


def apply_overview(args: argparse.Namespace, parser: argparse.ArgumentParser):
    """Expand ``--overview N`` into the flags it stands for, adding to any given explicitly."""
    n = args.overview
    if n < 1:
        parser.error("--overview needs a positive number, e.g. --overview 3")

    def extended(columns: list[str] | None, extra: str) -> list[str]:
        # A bare --highlight/--highlight-pic means the aspect ratios; keep them when adding
        if columns is None:
            return [extra]
        return (columns or ["ar1", "ar2"]) + [extra]

    args.highlight = extended(args.highlight, f"energies:{n}")
    args.highlight_pic = extended(args.highlight_pic, f"energies:{n}")
    args.colour_by = args.colour_by + [f"properties:{n}"]
    args.group = True

    LOG.info(
        "--overview %d: --group --highlight %s --highlight-pic %s --colour-by %s",
        n,
        " ".join(args.highlight),
        " ".join(args.highlight_pic),
        " ".join(args.colour_by),
    )
    if not args.energies and args.input_dir and not args.results_dir:
        LOG.warning("--overview: add --energies to read the interaction energies from net files")


def main(argv: list[str] | None = None) -> int:
    """Main entry point with command line interface"""
    parser = build_parser()
    args = parser.parse_args(argv)

    # Validate arguments
    if not any([args.input_dir, args.results_dir]):
        parser.error("Either --input_dir or --results_dir must be specified")

    if args.cda and not args.directions:
        parser.error("--directions must be specified for CDA analysis")

    if args.overview is not None:
        apply_overview(args, parser)

    # Set up configuration
    if args.output_dir:
        save_folder = args.output_dir
    elif args.input_dir:
        save_folder = args.input_dir / "RESULTS"
    else:
        save_folder = args.results_dir

    line_filter = {}
    for item in args.line_filter:
        col, _, raw_vals = item.partition("=")
        parsed = []
        for v in raw_vals.split(","):
            try:
                parsed.append(float(v))
            except ValueError:
                parsed.append(v)
        line_filter[col] = parsed

    zingg_mode = "bounding_box" if args.box else "svd"
    if args.box and not args.input_dir:
        LOG.warning(
            "--box has no effect when reading pre-computed results with --results_dir; "
            "Zingg ratios are only recalculated from structures in --input_dir",
        )
    # Shape files are recorded relative to the input folder. When re-plotting, the
    # results folder is by default <input_dir>/RESULTS, so its parent is the input folder.
    if args.input_dir:
        shapes_root = args.input_dir
    elif args.results_dir:
        shapes_root = args.results_dir.resolve().parent
    else:
        shapes_root = None

    config = ShapeAnalysisConfig(
        save_folder=save_folder,
        show_plots=args.show,
        lmax=args.lmax,
        zingg_method=zingg_mode,
        ar_limits=args.ar_limits,
        highlight=args.highlight,
        colour_by=args.colour_by,
        highlight_pictures=args.highlight_pic,
        input_dir=shapes_root,
        group_energies=args.group,
        zoom=args.zoom,
    )

    if args.group and not (args.solvent or args.cda):
        LOG.warning("--group only applies to energies read from net files (--solvent, --cda)")

    try:
        theme = build_theme(
            font_scale=args.font_scale,
            border_width=args.border_width,
            border_color=args.border_color,
        )
    except ValueError as exc:
        parser.error(str(exc))

    # Initialise pipeline
    pipeline = CrystalShapeAnalysisPipeline(config, plot_style=args.plot_style, custom_theme=theme)

    # Discover files if input directory provided
    files = {}
    if args.input_dir:
        files = FileDiscovery.find_files(args.input_dir)

    # Run analyses based on flags
    results = {}

    if not any([args.general, args.solvent, args.wulff, args.size, args.movies, args.cda]):
        LOG.error(
            "No analysis selected. Please provide at least one of: "
            "--general, --solvent, --wulff, --size, --movies, --cda",
        )
        return 1

    if args.general:
        if not args.input_dir:
            LOG.error("--input_dir required for general analysis")
            return 1

        # Use XYZ files for general analysis
        shape_files = files.get("xyz", [])
        if not shape_files:
            LOG.error("No XYZ files found for general analysis")
            return 1

        results["general"] = pipeline.run_general_analysis(
            shape_files, args.energy_csv[0] if args.energy_csv else None, args.labels
        )

    if args.solvent:
        if args.input_dir:
            # Use discovered files
            shape_files = files.get("xyz", [])
            occ_outputs = files.get("occ", []) if args.occ else None
        elif args.results_dir:
            # Load from existing CSV
            csv_path = args.results_dir / "cg_analysis.csv"
            if csv_path.exists():
                results["cg"] = pipeline.group_energies(pd.read_csv(csv_path))
                pipeline.plots.create_zingg_plot(results["cg"], name="cg")
                pipeline.plots.create_colored_zingg(results["cg"], name="cg", mode="parameter")
                pipeline.plots.create_annotated_zingg_plots(results["cg"], args.labels, name="cg")
                return 0

            LOG.error("CG results CSV not found at %s", csv_path)
            return 1

        if shape_files:
            results["solvent"] = pipeline.run_solvent_analysis(
                shape_files,
                args.solvent_json,
                occ_outputs,
                args.energies,
                args.exclude,
                args.labels,
            )

    if args.wulff:
        if not args.input_dir:
            LOG.error("--input_dir required for Wulff analysis")
            return 1

        # Use PLY files for Wulff analysis
        wulff_files = files.get("wulff", [])
        if not wulff_files:
            LOG.error("No PLY files found for Wulff analysis")
            return 1

        results["wulff"] = pipeline.run_wulff_analysis(
            wulff_files,
            args.energy_csv[0] if args.energy_csv else None,
            args.labels,
            line_filter=line_filter or None,
            smooth_window=args.smooth_window,
        )

    if args.size:
        if args.input_dir:
            size_files = files.get("size", [])
        elif args.results_dir:
            raise NotImplementedError(
                "Size analysis requires input path. Reading from results directory not supported yet."
            )

        if args.crystallography is None:
            raise FileNotFoundError(
                "Can't analyse size files without the crystallographic information! Provide JSON."
            )

        if not size_files:
            LOG.warning("No size files found to process")

        n_steps = args.wulff_interval

        try:
            all_results = process_multiple_size_files(
                size_files,
                Path(args.crystallography),
                save_folder,
                n_steps=n_steps,
                expand_symmetry=True,
                reduce_facets=True,
                save_data=True,
            )

            results["size"] = pipeline.run_size_wulff_analysis(
                all_results,
                args.energy_csv or None,
                relative=not args.no_relative,
                labels_to_show=args.labels,
                line_filter=line_filter or None,
                smooth_window=args.smooth_window,
                timecol=args.timecol,
                log_scale=args.log_scale,
            )

            LOG.info("Successfully processed %s size files", len(all_results))

        except Exception:
            LOG.exception("Error in Wulff shape generation workflow")

    if args.movies:
        if args.input_dir:
            movie_files = files.get("xyz", [])
        elif args.results_dir:
            raise NotImplementedError(
                "Size analysis requires input path. Reading from results directory not supported yet."
            )

        if not movie_files:
            LOG.warning("No movie (XYZ) files found to process")

        try:
            results["movies"] = pipeline.run_movie_analysis(
                movie_files,
                args.energy_csv[0] if args.energy_csv else None,
                args.labels,
                line_filter=line_filter or None,
                smooth_window=args.smooth_window,
            )

            LOG.info("Successfully processed %s movie files", len(results["movies"]))

        except Exception:
            LOG.exception("Error in Movies workflow")

    if args.cda:
        if args.input_dir:
            # Use discovered files
            cda_files = files.get("cda", [])
        elif args.results_dir:
            # Load from existing CSV
            csv_path = args.results_dir / "cda_analysis.csv"
            if csv_path.exists():
                results["cda"] = pipeline.group_energies(pd.read_csv(csv_path))
                pipeline.plots.create_zingg_plot(results["cda"], name="cda")
                pipeline.plots.create_annotated_zingg_plots(results["cda"], args.labels, name="cda")
                return 0

            LOG.error("CDA results CSV not found at %s", csv_path)
            return 1

        if cda_files and args.directions:
            results["cda"] = pipeline.run_cda_analysis(
                cda_files, args.directions, args.energies, args.labels
            )

    # Print summary
    if results:
        LOG.info("\n%s", "=" * 50)
        LOG.info("ANALYSIS SUMMARY")
        LOG.info("=" * 50)
        for analysis_type, df in results.items():
            LOG.info("%s: %d shapes analysed", analysis_type.upper(), len(df))
        LOG.info("\nResults saved to: %s", config.save_folder)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
