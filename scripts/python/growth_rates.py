#!/usr/bin/env python3
"""
Growth Rates

This script calculates and plots growth rates from size.csv files generated
by crystal growth simulations. It processes multiple simulation files,
calculates growth rates using linear fitting, and generates comprehensive
visualization plots.

Usage:
    python growth_rates.py --input /path/to/data --directions "100" "010" "001"

Requirements:
    - numpy
    - pandas
    - matplotlib
"""

import argparse
import logging
import re
from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("GrowthRatePlotter")


def get_x_axis(df: pd.DataFrame, *, time_col: str = "time", tol: float = 1e-12):
    """Return a suitable x-axis for growth-rate fitting.

    Uses the time column if it has a sufficient span; otherwise falls back
    to the row index.
    """
    if time_col not in df.columns:
        logger.debug("Using row indices as x-axis (time column not found)")
        return np.arange(len(df), dtype=float)

    x_time = df[time_col].to_numpy(dtype=float)

    if not np.isfinite(x_time).all() or len(x_time) < 2 or np.ptp(x_time) < tol:
        logger.debug("Using row indices as x-axis (time column unsuitable)")
        return np.arange(len(df), dtype=float)

    logger.debug("Using time column as x-axis")
    return x_time


def build_growthrates(
    size_file_list: List[Path],
    supersat_list: List[float],
    directions: List[str],
    plot_raw_data: bool = False,
    raw_data_output: Optional[Path] = None,
    xaxis_mode: str = "auto",
    time_tol: float = 1e-12,
) -> Optional[pd.DataFrame]:
    """
    Generate the growth rate dataframe from the size.csv files.

    Args:
        size_file_list: List of paths to size.csv files
        supersat_list: List of supersaturation values corresponding to each file
        directions: List of crystal direction strings (e.g., [" 1 0 0", " 0 1 0"])
        plot_raw_data: If True, plot raw size vs time data for each simulation
        raw_data_output: Directory to save raw data plots
        xaxis_mode: How to choose the x-axis for linear fitting:
            - "auto"  – use time column when valid, else fall back to row index for all files
            - "time"  – always use time column; files without a valid one are skipped
            - "index" – always use row index, ignoring any time column
        time_tol: Minimum time span to consider the time column valid

    Returns:
        DataFrame with columns: Simulation Number, Supersaturation, and growth rates per direction
    """
    n_size_files = len(size_file_list)

    if n_size_files == 0:
        logger.error("No size files provided")
        return None

    logger.info(f"{n_size_files} size files used to calculate growth rate data")
    logger.info(f"X-axis mode: {xaxis_mode}")

    if plot_raw_data and raw_data_output:
        raw_data_output = Path(raw_data_output)
        raw_data_output.mkdir(parents=True, exist_ok=True)
        logger.info(f"Raw data plots will be saved to {raw_data_output}")

    growth_list = []
    kept_supersats = []
    use_index_for_all = xaxis_mode == "index"
    restart = True

    while restart:
        restart = False
        growth_list = []
        kept_supersats = []

        for i, f in enumerate(size_file_list):
            f = Path(f)
            logger.info(f"Processing file {i + 1}/{n_size_files}: {f.name}")

            try:
                lt_df = pd.read_csv(f, encoding="utf-8", encoding_errors="replace")
            except Exception as e:
                logger.warning(f"Failed to read file {f.name}: {e}")
                continue

            # Check if all required directions are present
            missing_directions = [d for d in directions if d not in lt_df.columns]
            if missing_directions:
                logger.warning(
                    f"Skipping file {f.name}: missing direction columns {missing_directions}"
                )
                continue

            # Determine x-axis data
            if use_index_for_all:
                x_data = np.arange(len(lt_df), dtype=float)
            else:
                time_col = "time"
                if time_col not in lt_df.columns:
                    if xaxis_mode == "time":
                        logger.warning(
                            f"Skipping file {f.name}: time column not found (forced time mode)"
                        )
                        continue
                    # auto mode: restart using index for all files
                    logger.info(
                        f"Time column missing in file {f.name} - restarting with index for all files"
                    )
                    use_index_for_all = True
                    restart = True
                    break
                else:
                    x_time = lt_df[time_col].to_numpy(dtype=float)
                    if (
                        not np.isfinite(x_time).all()
                        or len(x_time) < 2
                        or np.ptp(x_time) < time_tol
                    ):
                        if xaxis_mode == "time":
                            logger.warning(
                                f"Skipping file {f.name}: time column unsuitable (forced time mode)"
                            )
                            continue
                        # auto mode: restart using index for all files
                        logger.info(
                            f"Time column unsuitable in file {f.name} - restarting with index for all files"
                        )
                        use_index_for_all = True
                        restart = True
                        break
                    else:
                        x_data = x_time

            tokens = re.findall(r"\d+", f.name)
            sim_num = int(tokens[-1]) if tokens else i

            # Keep rows only up to the first row where any direction is 0
            all_positive = np.all(
                [np.asarray(lt_df[d], dtype=float) > 0 for d in directions],
                axis=0,
            )
            first_false = np.argmin(all_positive)
            cutoff = first_false if not all_positive[first_false] else len(all_positive)
            mask = np.zeros(len(all_positive), dtype=bool)
            mask[:cutoff] = True

            # Plot raw data if requested
            if plot_raw_data and raw_data_output:
                fig, axes = plt.subplots(len(directions), 1, figsize=(10, 3 * len(directions)))
                if len(directions) == 1:
                    axes = [axes]

                supersat_val = supersat_list[i] if i < len(supersat_list) else "Unknown"
                fig.suptitle(
                    f"Sim {sim_num} - Supersaturation: {supersat_val} kcal/mol", fontsize=14
                )

                for idx, direction in enumerate(directions):
                    y_data = np.asarray(lt_df[direction], dtype=float)
                    x_plot = x_data[mask]
                    y_plot = y_data[mask]

                    if mask.sum() >= 2:
                        model = np.polyfit(x_plot, y_plot, 1)
                        fit_line = model[0] * x_plot + model[1]
                        axes[idx].plot(
                            x_plot, fit_line, "r-", linewidth=2, label=f"Fit (slope={model[0]:.6f})"
                        )

                    axes[idx].scatter(x_data, y_data, s=2, alpha=0.6, label="Data")
                    axes[idx].set_xlabel("Time" if not use_index_for_all else "Row Index")
                    axes[idx].set_ylabel(f"Size [{direction}]")
                    axes[idx].legend()
                    axes[idx].grid(True, alpha=0.3)
                    axes[idx].set_title(f"Direction: {direction}")

                plt.tight_layout()
                plot_filename = raw_data_output / f"raw_data_sim_{sim_num:03d}.png"
                plt.savefig(plot_filename, dpi=150)
                plt.close()
                logger.debug(f"Saved raw data plot: {plot_filename}")

            gr_list = [sim_num]
            for direction in directions:
                y_data = np.asarray(lt_df[direction], dtype=float)
                if mask.sum() < 2:
                    gr_list.append(0.0)
                    continue
                model = np.polyfit(x_data[mask], y_data[mask], 1)
                gr_list.append(model[0])

            growth_list.append(gr_list)
            kept_supersats.append(supersat_list[i])
            logger.info(f"Calculated growth rates for simulation {sim_num}")

    if not growth_list:
        logger.error("No valid data was processed")
        return None

    logger.debug(f"Growth Rate data (list): {growth_list}")
    growth_array = np.asarray(growth_list)
    gr_df = pd.DataFrame(growth_array, columns=["Simulation Number"] + directions)
    gr_df.insert(1, "Supersaturation", kept_supersats)

    return gr_df


def plot_growth_rates(gr_df: pd.DataFrame, directions: List[str], savepath: Path):
    """
    Generate comprehensive growth rate plots.

    Creates 6 different plots:
    1. Combined growth/dissolution rates
    2. Growth rates (supersaturation >= 0)
    3. Growth rates (zoomed view)
    4. Dissolution rates (supersaturation <= 0)
    5. Dissolution rates (zoomed view)

    Args:
        gr_df: DataFrame with growth rate data
        directions: List of direction column names to plot
        savepath: Directory path where plots will be saved
    """
    savepath = Path(savepath)
    savepath.mkdir(parents=True, exist_ok=True)

    # Sort by supersaturation for proper line plot connections
    gr_df = gr_df.sort_values("Supersaturation").reset_index(drop=True)
    logger.debug("Data sorted by supersaturation")

    x_data = gr_df["Supersaturation"]

    # Plot 1: Combined growth and dissolution rates
    plt.figure(figsize=(7, 5))
    for direction in directions:
        plt.scatter(x_data, gr_df[direction], s=1.2)
        plt.plot(x_data, gr_df[direction], label=direction)
    plt.legend()
    plt.xlabel("Supersaturation (kcal/mol)")
    plt.ylabel("Growth Rate")
    plt.tight_layout()
    logger.info("Plotting growth/dissolution rates")
    plt.savefig(savepath / "growth_and_dissolution_rates.png", dpi=300)
    plt.close()

    # Plot 2: Growth rates only (supersaturation >= 0)
    growth_data = gr_df[gr_df["Supersaturation"] >= 0].sort_values("Supersaturation")
    plt.figure(figsize=(5, 5))
    for direction in directions:
        plt.scatter(growth_data["Supersaturation"], growth_data[direction], s=1.2)
        plt.plot(growth_data["Supersaturation"], growth_data[direction], label=direction)
    plt.legend()
    plt.xlabel("Supersaturation (kcal/mol)")
    plt.ylabel("Growth Rate")
    plt.tight_layout()
    logger.info("Plotting growth rates")
    plt.savefig(savepath / "growth_rates.png", dpi=300)
    plt.close()

    # Plot 3: Growth rates (zoomed)
    plt.figure(figsize=(5, 5))
    for direction in directions:
        plt.scatter(growth_data["Supersaturation"], growth_data[direction], s=1.2)
        plt.plot(growth_data["Supersaturation"], growth_data[direction], label=direction)
    plt.legend()
    plt.xlabel("Supersaturation (kcal/mol)")
    plt.ylabel("Growth Rate")
    plt.xlim(0.0, 2.5)
    plt.ylim(0.0, 0.4)
    plt.tight_layout()
    logger.info("Plotting growth rates (zoomed)")
    plt.savefig(savepath / "growth_rates_zoomed.png", dpi=300)
    plt.close()

    # Plot 4: Dissolution rates (supersaturation <= 0)
    dissolution_data = gr_df[gr_df["Supersaturation"] <= 0].sort_values("Supersaturation")
    if not dissolution_data.empty:
        plt.figure(figsize=(7, 5))
        for direction in directions:
            plt.scatter(
                dissolution_data["Supersaturation"],
                dissolution_data[direction],
                label=direction,
                s=1.2,
            )
        plt.legend()
        plt.xlabel("Supersaturation (kcal/mol)")
        plt.ylabel("Dissolution Rate")
        plt.tight_layout()
        logger.info("Plotting dissolution rates")
        plt.savefig(savepath / "dissolution_rates.png", dpi=300)
        plt.close()

        # Plot 5: Dissolution rates (zoomed)
        plt.figure(figsize=(5, 5))
        for direction in directions:
            plt.scatter(dissolution_data["Supersaturation"], dissolution_data[direction], s=1.2)
            plt.plot(
                dissolution_data["Supersaturation"], dissolution_data[direction], label=direction
            )
        plt.legend()
        plt.xlabel("Supersaturation (kcal/mol)")
        plt.ylabel("Growth Rate")
        plt.xlim(-2.5, 0.0)
        plt.ylim(-2.5, 0.0)
        plt.tight_layout()
        logger.info("Plotting dissolution rates (zoomed)")
        plt.savefig(savepath / "dissolution_rates_zoomed.png", dpi=300)
        plt.close()
    else:
        logger.info("No dissolution data available (all supersaturation >= 0)")

    logger.info(f"All plots saved to {savepath}")


def extract_sim_number(filepath: Path) -> int:
    """Extract simulation number from file path."""
    tokens = re.findall(r"\d+", filepath.name)
    return int(tokens[-1]) if tokens else 0


def find_size_files_and_supersats(input_folder: Path) -> tuple[List[Path], List[float]]:
    """
    Find all size.csv files and extract supersaturation from simulation_parameters.txt files.

    Args:
        input_folder: Path to search for files

    Returns:
        Tuple of (size_files, supersaturation_values), both sorted by simulation number
    """
    input_folder = Path(input_folder)

    # Find all size files
    size_files = list(input_folder.rglob("*size.csv"))
    logger.info(f"Found {len(size_files)} size files")

    # Create a mapping of simulation number to size file and supersaturation
    sim_data = {}

    for size_file in size_files:
        sim_num = extract_sim_number(size_file)

        # Look for simulation_parameters.txt in the same directory
        param_file = list(size_file.parent.glob("*simulation_parameters.txt"))[0]
        supersat = None

        if param_file.exists():
            try:
                with open(param_file, "r", encoding="utf-8", errors="replace") as f:
                    for line in f:
                        if line.startswith("Starting delta mu value (kcal/mol):"):
                            supersat = float(line.split()[-1])
                            break
            except Exception as e:
                logger.warning(f"Failed to read {param_file}: {e}")

        if supersat is None:
            logger.warning(f"No supersaturation found for {size_file.name}, using 0.0")
            supersat = 0.0

        sim_data[sim_num] = {"size_file": size_file, "supersat": supersat}

    # Sort by simulation number
    sorted_sim_nums = sorted(sim_data.keys())

    sorted_size_files = [sim_data[num]["size_file"] for num in sorted_sim_nums]
    sorted_supersats = [sim_data[num]["supersat"] for num in sorted_sim_nums]

    logger.info(
        f"Extracted {len(sorted_supersats)} supersaturation values from simulation_parameters.txt files"
    )
    logger.info(
        f"Supersaturation range: {min(sorted_supersats):.2f} to {max(sorted_supersats):.2f} kcal/mol"
    )

    return sorted_size_files, sorted_supersats


def main():
    parser = argparse.ArgumentParser(
        description="Calculate and plot growth rates from crystal growth simulation data"
    )
    parser.add_argument(
        "--input", "-i", type=str, required=True, help="Input folder containing size.csv files"
    )
    parser.add_argument(
        "--directions",
        "-d",
        nargs="+",
        required=True,
        help='Crystal directions to analyze (e.g., " 1 0 0" " 0 1 0" " 0 0 1")',
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Output folder for plots and CSV (default: input_folder/growth_rate_analysis)",
    )
    parser.add_argument(
        "--supersats",
        "-s",
        nargs="+",
        type=float,
        default=None,
        help="Supersaturation values for each file (optional, will auto-detect if not provided)",
    )
    parser.add_argument(
        "--xaxis-mode",
        choices=["auto", "time", "index"],
        default="auto",
        help=(
            "X-axis mode for linear fitting: "
            "'auto' uses time when valid and falls back to row index, "
            "'time' forces time column (skips files without one), "
            "'index' always uses row index (default: auto)"
        ),
    )
    parser.add_argument(
        "--plot-raw",
        action="store_true",
        help="Plot raw size vs time data for each simulation to investigate fitting",
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Setup paths
    input_folder = Path(args.input)
    if not input_folder.exists():
        logger.error(f"Input folder does not exist: {input_folder}")
        return

    if args.output:
        output_folder = Path(args.output)
    else:
        output_folder = input_folder / "growth_rate_analysis"

    output_folder.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output will be saved to: {output_folder}")

    # Find size files and extract supersaturation values
    if args.supersats:
        # Manual supersaturation override - still need to find and sort size files
        logger.info("Using manually provided supersaturation values")
        size_files, _ = find_size_files_and_supersats(input_folder)
        supersat_list = args.supersats
        if len(supersat_list) != len(size_files):
            logger.warning(
                f"Number of supersaturation values ({len(supersat_list)}) "
                f"does not match number of files ({len(size_files)})"
            )
    else:
        # Automatic extraction from simulation_parameters.txt files
        logger.info("Auto-extracting supersaturation from simulation_parameters.txt files")
        size_files, supersat_list = find_size_files_and_supersats(input_folder)

    if not size_files:
        logger.error("No size files found!")
        return

    # Build growth rates dataframe
    logger.info(f"Calculating growth rates for directions: {args.directions}")

    # Setup raw data output folder if needed
    raw_data_folder = None
    if args.plot_raw:
        raw_data_folder = output_folder / "raw_data_plots"
        logger.info("Raw data plotting enabled")

    gr_df = build_growthrates(
        size_files,
        supersat_list,
        args.directions,
        plot_raw_data=args.plot_raw,
        raw_data_output=raw_data_folder,
        xaxis_mode=args.xaxis_mode,
    )

    if gr_df is None or gr_df.empty:
        logger.error("Failed to calculate growth rates")
        return

    # Save growth rates CSV
    csv_path = output_folder / "growthrates.csv"
    gr_df.to_csv(csv_path, index=False)
    logger.info(f"Growth rates saved to: {csv_path}")

    # Display DataFrame
    print("\nGrowth Rates DataFrame:")
    print(gr_df.to_string())

    # Generate plots
    logger.info("Generating growth rate summary plots...")
    plot_growth_rates(gr_df, args.directions, output_folder)

    if args.plot_raw:
        logger.info(f"Raw data plots saved to: {raw_data_folder}")

    logger.info("Done! Check the output folder for results.")


if __name__ == "__main__":
    main()
