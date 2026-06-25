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
from collections import defaultdict
from pathlib import Path
from typing import List, Optional, Dict, Tuple, Literal

import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import mplcursors
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from tqdm import tqdm
from natsort import natsorted

from cgpytools.shape_analysis import CrystalShape, ShapeAnalyser
from cgpytools.cg_net import CGNet
from cgpytools.plot import PlotTheme, GlobalPlotStyler
from cgpytools.surfaces import process_multiple_size_files
from cgpytools.log import setup_logging

LOG = setup_logging(name="CG-ANALYSE")


class ShapeAnalysisConfig:
    """Configuration class for shape analysis parameters"""

    def __init__(
        self,
        save_folder: Path,
        show_plots: bool = False,
        lmax: int = 20,
        zingg_method: Literal["bounding_box", "svd"] = "svd",
        ar_limits: bool = False,
    ):
        self.save_folder = save_folder
        self.show_plots = show_plots
        self.lmax = lmax
        self.zingg_method = zingg_method
        self.ar_limits = ar_limits

        # Create save folder if it doesn't exist
        self.save_folder.mkdir(parents=True, exist_ok=True)


class FileDiscovery:
    """Handles discovery and organization of input files"""

    @staticmethod
    def find_files(
        folderpath: Path,
        include_patterns: Dict[str, List[str]] | None = None,
        exclude_patterns: Dict[str, List[str]] | None = None,
    ) -> Dict[str, List[Path]]:
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

        file_types: Dict[str, List[Path]] = {}

        for category, patterns in include_patterns.items():
            found: List[Path] = []
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
    def parse_shape_name(value: str | int | float) -> str | int:
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
    ) -> Dict[str, List[Tuple[str, float]]]:
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

            LOG.info(f"Loaded energy data for {len(energy_data)} shapes from {csv_path}")
            return energy_data

        except Exception as e:
            LOG.error(f"Failed to load energy data from {csv_path}: {e}", exc_info=True)
            return {}

    @staticmethod
    def load_from_csvs(
        csv_paths: List[Path], shape_name_column: str = "shape_name"
    ) -> Dict[str, List[Tuple[str, float]]]:
        """Load and merge energy data from multiple CSV files"""
        merged: Dict[str, List[Tuple[str, float]]] = {}
        for p in csv_paths or []:
            if p and p.exists():
                merged.update(EnergyDataLoader.load_from_csv(p, shape_name_column))
        return merged

    @staticmethod
    def load_from_net_files(shapes: List[Path]) -> Dict[str, List[float]]:
        """Load energy data from net.txt files"""
        energy_data = {}

        for shape in shapes:
            shape = Path(shape)
            if shape.suffix == ".stdout":
                netfile = shape.parent / f"solvent_{shape.parent.name.split('_')[-1]}" / "net.txt"
            else:
                netfile = shape.parent / "net.txt"

            if netfile.exists():
                try:
                    net = CGNet(netfile)
                    net.parse()
                    energies = net.unique_energies_arr.flatten()
                    energies = energies[~np.isnan(energies)]  # Remove NaN values
                    energy_data[shape.stem] = energies.tolist()
                    LOG.debug(f"Loaded {len(energies)} energies from {netfile}")
                except Exception as e:
                    LOG.error(f"Failed to load energies from {netfile}: {e}")

        return energy_data


class SolventDataLoader:
    """Handles loading and processing of solvent-related data"""

    @staticmethod
    def load_solvent_properties(solvent_json: Path) -> Dict:
        """Load solvent properties from JSON file"""

        try:
            if not Path(solvent_json).is_file or Path(solvent_json).stat().st_size == 0:
                LOG.error(f"File is empty/not found: {solvent_json}")
                return FileNotFoundError
            with open(solvent_json, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            LOG.error(f"Failed to load solvent properties from {solvent_json}: {e}", exc_info=True)
            return {}

    @staticmethod
    def extract_solvent_name(file_path: Path) -> Optional[str]:
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
    def load_occ_solubility(occ_outputs: List[Path]) -> Dict[str, List[float]]:
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
                    LOG.info(f"Solubility in {solvent}: {solubility}")
                    occ_info[solvent].append(solubility)
                else:
                    LOG.warning(f"Solubility not found in {output}")

            except Exception as e:
                LOG.error(f"Failed to process OCC output {output}: {e}")

        return dict(occ_info)


class CrystalShapeAnalyser:
    """Main class for analyzing crystal shapes"""

    def __init__(self, config: ShapeAnalysisConfig):
        self.config = config
        self.shape_analyser = ShapeAnalyser(zingg_method=self.config.zingg_method)

    def analyse_general_shapes(
        self,
        shape_files: List[Path],
        energy_csv: Optional[Path] = None,
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
                shape_info["ar1"].append(metrics.aspect1)
                shape_info["ar2"].append(metrics.aspect2)
                shape_info["sa"].append(metrics.surface_area)
                shape_info["vol"].append(metrics.volume)
                shape_info["sa_vol"].append(metrics.surface_area_to_volume_ratio)

                if shape_name in energy_data:
                    for i, (col, energy) in enumerate(energy_data[shape_name], 1):
                        shape_info[col].append(energy)

                LOG.debug(f"Analysed shape: {shape_name}")

            except Exception as e:
                LOG.error(f"Failed to analyse shape {shape_file}: {e}")
                continue
        LOG.info("\n")

        return self._finalise_dataframe(shape_info, "general_shapes")

    def analyse_movies(
        self,
        movie_files: List[Path],
        energy_csv: Optional[Path] = None,
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

                for frame_idx, metrics in self.shape_analyser.get_all_frame_metrics().items():
                    movie_info["shape_name"].append(movie_name)
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

                LOG.debug(f"Analysed movie: {movie_name} ({len(crystal.frames)} frames)")

            except Exception as e:
                LOG.error(f"Failed to analyse movie {movie_file}: {e}")
                continue
        LOG.info("\n")

        return self._finalise_dataframe(movie_info, "movies")

    def analyse_solvent_shapes(
        self,
        shapes: List[Path],
        solvent_json: Path,
        occ_outputs: Optional[List[Path]] = None,
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
                LOG.error(f"Couldn't find solvent for {shape_file}")
                continue
            try:
                crystal = CrystalShape.from_file(shape_file)
                self.shape_analyser.analyse_crystal(
                    crystal,
                    frame_idx=-1,
                )
                metrics = self.shape_analyser.get_frame_metrics(-1)
                shape_info["solvent"].append(solvent)
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
            except Exception as e:
                LOG.error(f"Failed to analyse shape {shape_file}: {e}", exc_info=True)
                continue

        LOG.info("\n")

        return self._finalise_dataframe(shape_info, "solvent_shapes")

    def analyse_cda_shapes(
        self, cda_files: List[Path], directions: List[str], get_energies: bool = False
    ) -> pd.DataFrame:
        """Analyse CDA simulation results"""
        ar_dict = defaultdict(list)

        for cda_file in tqdm(cda_files, desc="Analyzing CDA files"):
            solvent = SolventDataLoader.extract_solvent_name(cda_file)
            if not solvent:
                LOG.warning(f"Could not extract solvent name from {cda_file}")
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
                    LOG.warning(f"Frame information not found in {cda_file}")
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
                        except Exception as e:
                            LOG.error(f"Failed to load energies from {netfile}: {e}")

            except Exception as e:
                LOG.error(f"Failed to process CDA file {cda_file}: {e}")
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
        generation_results: Dict[str, Dict],
        energy_csvs: Optional[List[Path]] = None,
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

                    LOG.debug(f"Analyzed Wulff shape: {shape_file.name}")

                except Exception as e:
                    LOG.error(f"Failed to analyze shape {size_file} at row {i}: {e}")
                    continue

        if not shape_info:
            LOG.warning("No shapes were successfully analyzed")
            return pd.DataFrame()

        # Convert to DataFrame and return
        df = pd.DataFrame(shape_info)
        LOG.info(f"Successfully analyzed {len(df)} Wulff shapes")
        df.to_csv(self.config.save_folder / "size_analysis.csv", index=False)

        return df

    def _get_solvent_from_shape(self, shape_file: Path) -> Optional[str]:
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
            except Exception as e:
                LOG.error(f"Failed to read solvent from {shape_file}: {e}")

        return solvent

    def _finalise_dataframe(self, shape_info: Dict, name: str) -> pd.DataFrame:
        """Finalise the shape information dictionary into a DataFrame"""
        # Find the maximum list length
        max_length = max(len(lst) for lst in shape_info.values() if lst)

        # Pad shorter lists with zeros and warn
        for key, value in shape_info.items():
            if len(value) < max_length:
                original_length = len(value)
                shape_info[key] = value + [0] * (max_length - len(value))
                LOG.warning(f"{key} values were padded with {max_length - original_length} zeros")

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
    def get_subplot_layout(n: int) -> Tuple[int, int]:
        """Calculate optimal subplot layout"""
        n_cols = int(math.ceil(math.sqrt(n)))
        n_rows = int(math.ceil(n / n_cols))
        return n_rows, n_cols

    def _apply_ar_limits(self, ax, x_is_ar=True, y_is_ar=True):
        """Apply aspect ratio limits if configured"""
        if self.config.ar_limits:
            if x_is_ar:
                ax.set_xlim(0, 1)
            if y_is_ar:
                ax.set_ylim(0, 1)

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

        LOG.info(f"Plotting Basic Zingg Plot ({name})")

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
            x_cols = [col for col in df_clean.columns if col.startswith("Int_")]
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
            int_cols = [col for col in df_clean.columns if col.startswith("Int_")]
            exclude = set(common_cols + int_cols)
            x_cols = [col for col in df_clean.columns if col not in exclude]
        else:
            LOG.error(f"Unknown mode '{mode}'. Use 'energy' or 'parameter'.")
            return

        if not x_cols:
            LOG.warning(f"No x columns found for mode={mode}")
            return

        # Check that all z_vars exist in the dataframe
        missing_z_vars = [z for z in z_vars if z not in df_clean.columns]
        if missing_z_vars:
            LOG.error(f"Missing z variables in dataframe: {missing_z_vars}")
            return

        # Create separate figures for each z variable
        has_subfolders = "subfolder" in df_clean.columns
        for z_var in z_vars:
            LOG.info(f"Plotting Heatmap Plot for {mode.title()} Variables vs {z_var} ({name})")
            if has_subfolders:
                unique_subfolders = sorted(df_clean["subfolder"].dropna().unique())
                for x_var in x_cols:
                    self._create_subfolder_heatmap_figure(
                        df_clean, x_var, unique_subfolders, z_var, name, timecol, smooth_window=smooth_window, log_scale=log_scale
                    )
            else:
                self._create_heatmap_figure(df_clean, x_cols, z_var, name, mode, timecol, smooth_window=smooth_window, log_scale=log_scale)

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
        x_label = r"$\mu$" if x_var.lower() in ("mu", "supersat", "supersaturation") else x_var.upper()

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
                    fig.colorbar(hb, ax=axes[i], label=f"sign\u00b7log\u2081\u2080(|{z_display}|)" if log_scale else z_display, shrink=0.8)
                else:
                    pivot_data = subset.groupby([x_var, timecol])[z_var].mean().reset_index()
                    heatmap_data = pivot_data.pivot(index=x_var, columns=timecol, values=z_var)
                    if log_scale:
                        heatmap_data = np.sign(heatmap_data) * np.log10(np.abs(heatmap_data).clip(lower=1e-30))

                    vmin = None
                    vmax = None
                    heatmap_center = heatmap_data.mean().mean()
                    if self.config.ar_limits and z_var in ["ar1", "ar2"]:
                        vmin = 0
                        vmax = 1
                        heatmap_center = None

                    log_label = f"sign\u00b7log\u2081\u2080(|{z_display}|)" if log_scale else z_display
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

            except Exception as e:
                LOG.warning(f"Could not create heatmap for subfolder {subfolder}: {e}")
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
                    fig.colorbar(hb, ax=axes[i], label=f"sign\u00b7log\u2081\u2080(|{z_display}|)" if log_scale else z_display, shrink=0.8)
                else:
                    pivot_data = df_clean.groupby([x_var, timecol])[z_var].mean().reset_index()
                    heatmap_data = pivot_data.pivot(index=x_var, columns=timecol, values=z_var)
                    if log_scale:
                        heatmap_data = np.sign(heatmap_data) * np.log10(np.abs(heatmap_data).clip(lower=1e-30))

                    vmin = None
                    vmax = None
                    heatmap_center = heatmap_data.mean().mean()
                    if self.config.ar_limits and z_var in ["ar1", "ar2"]:
                        vmin = 0
                        vmax = 1
                        heatmap_center = None

                    log_label = f"sign\u00b7log\u2081\u2080(|{z_display}|)" if log_scale else z_display
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
                axes[i].set_ylabel(r"$\mu$" if x_var.lower() in ("mu", "supersat", "supersaturation") else x_var.upper())

            except Exception as e:
                LOG.warning(f"Could not create heatmap for {x_var}: {e}")
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
        line_filter: Optional[Dict[str, List]] = None,
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
            x_cols = [col for col in df_clean.columns if col.startswith("Int_")]
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
            int_cols = [col for col in df_clean.columns if col.startswith("Int_")]
            exclude = set(common_cols + int_cols)
            x_cols = [col for col in df_clean.columns if col not in exclude]
        else:
            LOG.error(f"Unknown mode '{mode}'. Use 'energy' or 'parameter'.")
            return

        if not x_cols:
            LOG.warning(f"No x columns found for mode={mode}")
            return

        missing_z_vars = [z for z in z_vars if z not in df_clean.columns]
        if missing_z_vars:
            LOG.error(f"Missing z variables in dataframe: {missing_z_vars}")
            return

        has_subfolders = "subfolder" in df_clean.columns
        for z_var in z_vars:
            LOG.info(f"Plotting Line Plot for {mode.title()} Variables vs {z_var} ({name})")
            if has_subfolders:
                unique_subfolders = sorted(df_clean["subfolder"].dropna().unique())
                for x_var in x_cols:
                    self._create_subfolder_lineplot_figure(
                        df_clean, x_var, unique_subfolders, z_var, name, mode, timecol, line_filter, smooth_window, log_scale=log_scale
                    )
            else:
                self._create_lineplot_figure(
                    df_clean, x_cols, z_var, name, mode, timecol, line_filter, smooth_window, log_scale=log_scale
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
        line_filter: Optional[Dict[str, List]] = None,
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
        x_label = r"$\mu$" if x_var.lower() in ("mu", "supersat", "supersaturation") else x_var.upper()

        # Collect all unique variable values across all subfolders
        all_x_vals = sorted(df_clean[x_var].dropna().unique())
        if line_filter and x_var in line_filter:
            allowed = set(line_filter[x_var])
            all_x_vals = [v for v in all_x_vals if v in allowed]

        if not all_x_vals:
            LOG.warning(f"No values to plot for {x_var}")
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
        line_filter: Optional[Dict[str, List]] = None,
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
                x_label = r"$\mu$" if x_var.lower() in ("mu", "supersat", "supersaturation") else x_var.upper()
                use_colorbar = not (line_filter and x_var in line_filter)

                if use_colorbar:
                    all_numeric = all(isinstance(v, (int, float)) for v in unique_x_vals)
                    if all_numeric:
                        val_arr = np.array(unique_x_vals, dtype=float)
                        norm = mcolors.Normalize(vmin=val_arr.min(), vmax=val_arr.max())
                        cmap = matplotlib.colormaps["viridis"]

                        def color_for(val):
                            return cmap(norm(float(val)))
                    else:
                        cmap = matplotlib.colormaps["tab10"]
                        val_index = {v: idx for idx, v in enumerate(unique_x_vals)}
                        norm = mcolors.BoundaryNorm(
                            boundaries=np.arange(-0.5, len(unique_x_vals) + 0.5),
                            ncolors=len(unique_x_vals),
                        )

                        def color_for(val):
                            return cmap(val_index[val] / max(len(unique_x_vals) - 1, 1))
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

            except Exception as e:
                LOG.warning(f"Could not create line plot for {x_var}: {e}")
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
            LOG.error(f"Missing required columns: {missing_cols}")
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
            LOG.warning(f"Color variable '{color_var}' not found in {mode} mode. Using timestep.")
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
            cbar_label = "Energy (kJ/mol)" if color_col.startswith("Int_") else color_col
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
                colorbar_title = "Energy (kJ/mol)" if color_col.startswith("Int_") else color_col
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
                    marker=dict(
                        size=8,
                        color=df_clean[color_col],
                        colorscale="Viridis",
                        colorbar=dict(
                            title=dict(text=colorbar_title, side="right"),
                            x=1.05,  # push right to avoid overlap with legend
                            y=0.5,  # vertical center
                            len=0.75,  # shrink length
                        ),
                        line=dict(width=1, color="black"),
                    ),
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
                    line=dict(color="red", width=4, dash="dash"),
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
                    line=dict(color="red", width=4, dash="dash"),
                    name="Zingg Boundary (Horizontal)",
                    showlegend=True,
                )
            )

            # Update layout
            x_label = "S:M" if x_col == "ar1" else x_col
            y_label = "M:L" if y_col == "ar2" else y_col

            # Build scene dict
            scene_dict = dict(
                xaxis_title=x_label,
                yaxis_title=y_label,
                zaxis_title=timecol.title(),
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.2)),
            )

            if self.config.ar_limits:
                scene_dict["xaxis"] = dict(title=x_label, range=[0, 1])
                scene_dict["yaxis"] = dict(title=y_label, range=[0, 1])

            fig.update_layout(
                title=f"Interactive 3D Zingg Evolution - {name.upper()} ({mode.title()} Mode)",
                scene=scene_dict,
                width=1000,
                height=800,
                legend=dict(
                    x=0.02,
                    y=0.98,  # top-left corner
                    bgcolor="rgba(255,255,255,0.6)",
                    bordercolor="black",
                    borderwidth=1,
                ),
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
            color_cols = [col for col in df_clean.columns if col.startswith("Int_")]
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
            int_cols = [col for col in df_clean.columns if col.startswith("Int_")]
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
            LOG.error(f"Unknown mode '{mode}'. Use 'energy' or 'parameter'.")
            return

        # Create plots for each color variable
        for i, color_var in enumerate(color_vars):
            LOG.info(f"Creating 3D Zingg plot colored by {color_var} ({mode} mode)")

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
            color_cols = [col for col in df_clean.columns if col.startswith("Int_")]
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
            int_cols = [col for col in df_clean.columns if col.startswith("Int_")]
            exclude = set(common_cols + int_cols)
            color_cols = [col for col in df_clean.columns if col not in exclude]
            cbar_label = "Parameter Value"
        else:
            LOG.error(f"Unknown mode '{mode}'. Use 'energy' or 'parameter'.")
            return

        if not color_cols:
            LOG.warning(f"No columns found for mode={mode}")
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

        LOG.info(f"Plotting {mode.title()} Coloured Zingg Plots ({name})")

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

        LOG.info(f"Plotting Correlation Matrices ({name})")

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

    def plot_solvent_effects(self, df: pd.DataFrame, name: str = "", exclude: List[str] = None):
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
        common_cols = ["solvent", "ar1", "ar2", "sa", "vol", "sa_vol", "shape_name"]
        int_cols = [col for col in df_clean.columns if col.startswith("Int_")]
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

        LOG.info(f"Plotting Solvent Parameter Effects on Shape({name})")

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
        energy_cols = [col for col in df_clean.columns if col.startswith("Int_")]

        if not energy_cols:
            LOG.warning("No interaction energy columns found")
            return

        LOG.info(f"Plotting Solvent Parameter Effects on Energies ({name})")

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

    def create_labeled_zingg_plot(
        self, df: pd.DataFrame, labels_to_show: List[str], name: str = ""
    ):
        """Create Zingg plot with specific labels highlighted using consistent styling"""
        df_clean = self._clean_columns(df)

        fig, ax = self.styler.create_figure(figsize=(12, 10))

        # Determine aspect ratio columns
        if name == "cda":
            ar_cols = [col for col in df_clean.columns if col.startswith("AspectRatio")]
            x, y = df_clean[ar_cols[0]], df_clean[ar_cols[1]]
            xlabel, ylabel = ar_cols[0], ar_cols[1]
        else:
            x, y = df_clean["ar1"], df_clean["ar2"]
            xlabel, ylabel = "S:M", "M:L"

        LOG.info(f"Plotting Labelled Zingg Diagrams ({name})")
        LOG.info(f"Labels: {labels_to_show}")

        # Main scatter plot
        colors = self.styler.get_color_palette(1)[0]
        ax.scatter(
            x,
            y,
            alpha=self.styler.theme.alpha_scatter,
            s=self.styler.theme.marker_size,
            color=colors,
        )

        # Apply Zingg styling
        title = f"Labelled Zingg Plot - {name.upper()}" if name else "Labelled Zingg Plot"
        self.styler.apply_zingg_style(ax, title=title, show_legend=False)

        # Override labels for CDA plots
        if name == "cda":
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

        # Add labels for specified items
        label_column = "solvent" if "solvent" in df_clean.columns else "shape_name"
        n_labels = len(labels_to_show)
        half = n_labels // 2

        # Get contrasting colors for labels
        self.styler.get_color_palette(2)
        highlight_color = self.styler.theme.primary_colors["orange"]

        for idx, row in df_clean.iterrows():
            label_value = row[label_column]
            if label_value in labels_to_show:
                if label_value in labels_to_show[:half]:
                    # Position labels above and to the right
                    xytext = (100, 100)
                    connection_style = "arc,angleA=90,angleB=0,armA=0,armB=-100,rad=0"
                else:
                    # Position labels below and to the right
                    xytext = (150, 5)
                    connection_style = "arc,angleA=-90,angleB=0,armA=0,armB=-160,rad=0"

                ax.annotate(
                    label_value,
                    (x.iloc[idx], y.iloc[idx]),
                    xytext=xytext,
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    bbox=dict(
                        boxstyle="round,pad=0.3",
                        fc=highlight_color,
                        alpha=0.8,
                        edgecolor="white",
                        linewidth=1,
                    ),
                    arrowprops=dict(
                        arrowstyle="->",
                        connectionstyle=connection_style,
                        color=self.styler.theme.primary_colors["red"],
                        alpha=0.8,
                        linewidth=1.5,
                    ),
                    fontsize=self.styler.theme.font_size_tick,
                    fontweight="medium",
                )

        plt.tight_layout()
        save_path = self.config.save_folder / f"zingg_labeled_{name}.png"
        self.styler.save_figure(fig, save_path)

        if self.config.show_plots:
            plt.show()
        else:
            plt.close()

    def set_plot_style(self, style: str, custom_theme: PlotTheme = None):
        """Change the global plot style"""
        self.styler.reset_style()
        self.styler = GlobalPlotStyler(theme=custom_theme, style=style)
        LOG.info(f"Plot style changed to: {style}")

    def __del__(self):
        """Cleanup: reset matplotlib settings when object is destroyed"""
        if hasattr(self, "styler"):
            self.styler.reset_style()


class CrystalShapeAnalysisPipeline:
    """Main pipeline class that orchestrates the entire analysis"""

    def __init__(self, config: ShapeAnalysisConfig):
        self.config = config
        self.crystal_analyser = CrystalShapeAnalyser(config)
        self.plots = ShapePlots(config)

    def run_general_analysis(
        self,
        shape_files: List[Path],
        energy_csv: Optional[Path] = None,
        labels_to_show: Optional[List[str]] = None,
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

        # Create labeled plot if labels specified
        if labels_to_show:
            self.plots.create_labeled_zingg_plot(df, labels_to_show, name="general")

        # If energy data is available, create energy-colored plots
        param_cols = [col for col in df.columns if col.startswith(("Int_", "x_", "X_"))]
        if param_cols:
            self.plots.create_colored_zingg(df, name="general", mode="parameter")

        LOG.info("General shape analysis completed")
        return df

    def run_solvent_analysis(
        self,
        shape_files: List[Path],
        solvent_json: Path,
        occ_outputs: Optional[List[Path]] = None,
        get_energies: bool = False,
        exclude_solvents: Optional[List[str]] = None,
        labels_to_show: Optional[List[str]] = None,
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

        # Create visualisations
        self.plots.create_zingg_plot(df, name="solvent")
        self.plots.create_colored_zingg(df, name="solvent", mode="parameter")
        self.plots.create_correlation_matrix(df, name="solvent")
        self.plots.plot_solvent_effects(df, name="solvent", exclude=exclude_solvents)

        # Create labeled plot if labels specified
        if labels_to_show:
            self.plots.create_labeled_zingg_plot(df, labels_to_show, name="solvent")

        # Create interaction energy plots if available
        energy_cols = [col for col in df.columns if col.startswith("Int_")]
        if energy_cols:
            self.plots.plot_interaction_energies(df, name="solvent")
            self.plots.create_colored_zingg(df, name="solvent", mode="energy")

        LOG.info("Solvent screening analysis completed")
        return df

    def run_cda_analysis(
        self,
        cda_files: List[Path],
        directions: List[str],
        get_energies: bool = False,
        labels_to_show: Optional[List[str]] = None,
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

        # Create visualisations
        self.plots.create_zingg_plot(df, name="cda")

        # Create labeled plot if labels specified
        if labels_to_show:
            self.plots.create_labeled_zingg_plot(df, labels_to_show, name="cda")

        # Create energy plots if available
        energy_cols = [col for col in df.columns if col.startswith(["Int_"])]
        if energy_cols:
            self.plots.create_colored_zingg(df, name="cda", mode="energy")

        LOG.info("CDA analysis completed")
        return df

    def run_size_wulff_analysis(
        self,
        generation_results: Dict[str, Dict],
        energy_csvs: Optional[List[Path]] = None,
        relative: bool = True,
        labels_to_show: Optional[List[str]] = None,
        line_filter: Optional[Dict[str, List]] = None,
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

        # Create labeled plot if labels specified
        if labels_to_show:
            self.plots.create_labeled_zingg_plot(df, labels_to_show, name="movies")

        self.plots.plot_var_heatmaps(df, name="wulff", mode="parameter", z_vars=["ar1", "ar2"], timecol=timecol, smooth_window=smooth_window)
        self.plots.plot_var_heatmaps(df, name="wulff", mode="parameter", z_vars=["sa", "vol"], timecol=timecol, smooth_window=smooth_window, log_scale=log_scale)
        self.plots.plot_var_heatmaps(df, name="wulff", mode="parameter", z_vars=["sa_vol"], timecol=timecol, smooth_window=smooth_window)
        self.plots.plot_var_lineplots(
            df, name="wulff", mode="parameter", z_vars=["ar1", "ar2"], timecol=timecol, line_filter=line_filter, smooth_window=smooth_window
        )
        self.plots.plot_var_lineplots(
            df, name="wulff", mode="parameter", z_vars=["sa", "vol"], timecol=timecol, line_filter=line_filter, smooth_window=smooth_window, log_scale=log_scale
        )
        self.plots.plot_var_lineplots(
            df, name="wulff", mode="parameter", z_vars=["sa_vol"], timecol=timecol, line_filter=line_filter, smooth_window=smooth_window
        )

        self.plots.create_multiple_3d_zingg_plots(df, name="wulff", mode="parameter")

        # Check if energy data is available

        energy_cols = [col for col in df.columns if col.startswith("Int_")]
        if energy_cols:
            # You could also do energy columns vs timestep with different z variables
            self.plots.create_colored_zingg(df, name="wulff", mode="energy")
            self.plots.plot_var_heatmaps(df, name="wulff", mode="energy", z_vars=["sa_vol"], timecol=timecol, smooth_window=smooth_window)
            self.plots.plot_var_lineplots(
                df, name="wulff", mode="energy", z_vars=["sa_vol"], timecol=timecol, line_filter=line_filter, smooth_window=smooth_window
            )
            self.plots.create_multiple_3d_zingg_plots(df, name="wulff", mode="energy")
        else:
            LOG.warning("No energy data available, skipping potential plots")

        LOG.info(f"Size/Wulff shape analysis completed for {len(df)} data points")
        return df

    def run_movie_analysis(
        self,
        movie_files: List[Path],
        energy_csv: Optional[Path] = None,
        labels_to_show: Optional[List[str]] = None,
        line_filter: Optional[Dict[str, List]] = None,
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

        # Create labeled plot if labels specified
        if labels_to_show:
            self.plots.create_labeled_zingg_plot(df, labels_to_show, name="movies")

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

        energy_cols = [col for col in df.columns if col.startswith("Int_")]
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
            f"Movie shape analysis completed for {len(df)} frames across {df['shape_name'].nunique()} movies"
        )
        return df

    def run_wulff_analysis(
        self,
        wulff_files: List[Path],
        energy_csv: Optional[Path] = None,
        labels_to_show: Optional[List[str]] = None,
        line_filter: Optional[Dict[str, List]] = None,
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

        # Create labeled plot if labels specified
        if labels_to_show:
            self.plots.create_labeled_zingg_plot(df, labels_to_show, name="wulff_ply")

        # If energy data is available, create energy-colored plots
        param_cols = [col for col in df.columns if col.startswith(("Int_", "x_", "X_"))]
        if param_cols:
            self.plots.create_colored_zingg(df, name="wulff_ply", mode="parameter")

        LOG.info("Wulff PLY shape analysis completed")
        return df


def main():
    """Main entry point with command line interface"""
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
        "-l", "--labels", nargs="+", default=[], help="Specific shapes/solvents to label in plots"
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
        "--timecol",
        choices=["timestep", "time", "frame"],
        default="timestep",
        help="Column to use as the time axis in line/heat plots (default: timestep). "
             "Use 'time' for the real time value from size.csv.",
    )

    args = parser.parse_args()

    # Validate arguments
    if not any([args.input_dir, args.results_dir]):
        raise ValueError("Either --input_dir or --results_dir must be specified")

    if args.cda and not args.directions:
        raise ValueError("--directions must be specified for CDA analysis")

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

    zingg_mode = "svd" if not args.box else "bounding_box"
    config = ShapeAnalysisConfig(
        save_folder=save_folder,
        show_plots=args.show,
        lmax=args.lmax,
        zingg_method=zingg_mode,
        ar_limits=args.ar_limits,
    )

    # Initialise pipeline
    pipeline = CrystalShapeAnalysisPipeline(config)

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

    if args.general:
        if not args.input_dir:
            LOG.error("--input_dir required for general analysis")
            return

        # Use XYZ files for general analysis
        shape_files = files.get("xyz", [])
        if not shape_files:
            LOG.error("No XYZ files found for general analysis")
            return

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
                results["cg"] = pd.read_csv(csv_path)
                pipeline.plots.create_zingg_plot(results["cg"], name="cg")
                pipeline.plots.create_colored_zingg(results["cg"], name="cg", mode="parameter")
                if args.labels:
                    pipeline.plots.create_labeled_zingg_plot(results["cg"], args.labels, name="cg")
            else:
                LOG.error(f"CG results CSV not found at {csv_path}")
            return

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
            return

        # Use PLY files for Wulff analysis
        wulff_files = files.get("wulff", [])
        if not wulff_files:
            LOG.error("No PLY files found for Wulff analysis")
            return

        results["wulff"] = pipeline.run_wulff_analysis(
            wulff_files, args.energy_csv[0] if args.energy_csv else None, args.labels,
            line_filter=line_filter or None, smooth_window=args.smooth_window,
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

            LOG.info(f"Successfully processed {len(all_results)} size files")

        except Exception as e:
            LOG.error(f"Error in Wulff shape generation workflow: {e}")

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

            LOG.info(f"Successfully processed {len(results['movies'])} movie files")

        except Exception as e:
            LOG.error(f"Error in Movies workflow: {e}")

    if args.cda:
        if args.input_dir:
            # Use discovered files
            cda_files = files.get("cda", [])
        elif args.results_dir:
            # Load from existing CSV
            csv_path = args.results_dir / "cda_analysis.csv"
            if csv_path.exists():
                results["cda"] = pd.read_csv(csv_path)
                pipeline.plots.create_zingg_plot(results["cda"], name="cda")
                if args.labels:
                    pipeline.plots.create_labeled_zingg_plot(
                        results["cda"], args.labels, name="cda"
                    )
            else:
                LOG.error(f"CDA results CSV not found at {csv_path}")
            return

        if cda_files and args.directions:
            results["cda"] = pipeline.run_cda_analysis(
                cda_files, args.directions, args.energies, args.labels
            )

    # Print summary
    if results:
        LOG.info("\n" + "=" * 50)
        LOG.info("ANALYSIS SUMMARY")
        LOG.info("=" * 50)
        for analysis_type, df in results.items():
            LOG.info(f"{analysis_type.upper()}: {len(df)} shapes analysed")
        LOG.info(f"\nResults saved to: {config.save_folder}")


if __name__ == "__main__":
    main()
