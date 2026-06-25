#!/usr/bin/env python3
"""
Growth Kinetics

Plots size_analysis.csv data from CG simulation subfolders.
For each subfolder, scatter-plots time (or timestep) vs a chosen parameter,
with a moving-average line per supersaturation value (Δμ), coloured by a
shared colorbar.

Optionally (--rates) opens a second figure with one panel per subfolder
showing Δμ (x) vs growth rate (y), where the growth rate is the linear
slope of the parameter over the x-axis for each supersaturation track.

Usage:
    python growth_kinetics.py --root /path/to/simulation --param vol
    python growth_kinetics.py --param ar1 --xaxis timestep --window 20
    python growth_kinetics.py --param vol --supersats -1.0 -2.0 -3.0
    python growth_kinetics.py --param vol --rates
    python growth_kinetics.py --param vol --rates --output fig.png

Requirements:
    numpy, pandas, matplotlib, scipy
"""

import argparse
import sys
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Columns present in size_analysis.csv
VALID_PARAMS = ["ar1", "ar2", "sa", "vol", "sa_vol"]
XAXIS_OPTIONS = ["time", "timestep"]

PARAM_LABELS = {
    "ar1": "Aspect Ratio 1",
    "ar2": "Aspect Ratio 2",
    "sa": "Surface Area",
    "vol": "Volume",
    "sa_vol": "SA / Volume",
}

XAXIS_LABELS = {
    "time": "Time",
    "timestep": "Timestep",
}


def sigmoid(x, A, k, x0, C):
    return A / (1 + np.exp(-k * (x - x0))) + C


def fit_sigmoid_x0(x, y):
    """
    Fit sigmoid and return x0 (transition point).
    Falls back safely if fit fails.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    # --- initial guesses (very important for stability) ---
    A0 = np.max(y) - np.min(y)
    C0 = np.min(y)
    x0_0 = x[np.argmin(np.abs(y - np.median(y)))]  # rough center
    k0 = 1.0 / (np.std(x) + 1e-8)

    p0 = [A0, k0, x0_0, C0]

    try:
        popt, _ = curve_fit(sigmoid, x, y, p0=p0, maxfev=10000)
        _, _, x0, _ = popt
        return x0
    except RuntimeError:
        # fallback: plateau-based estimate
        idx = np.argmin(np.abs(y))
        return x[idx]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(description="Plot MC time evolution from size_analysis.csv files.")
    p.add_argument(
        "--root",
        "-r",
        type=Path,
        default=Path("."),
        help=(
            "Root directory to search for subfolders containing "
            "RESULTS/size_analysis.csv (default: current directory)"
        ),
    )
    p.add_argument(
        "--param",
        "-p",
        default="vol",
        choices=VALID_PARAMS,
        help="Y-axis parameter to plot (default: vol)",
    )
    p.add_argument(
        "--xaxis",
        "-x",
        default="time",
        choices=XAXIS_OPTIONS,
        help="X-axis column: 'time' or 'timestep' (default: time)",
    )
    p.add_argument(
        "--window",
        "-w",
        type=int,
        default=None,
        metavar="N",
        help="Moving-average window size in number of points. If omitted, no line is drawn.",
    )
    p.add_argument(
        "--supersats",
        "-s",
        type=float,
        nargs="+",
        default=None,
        metavar="SUPERSAT",
        help=(
            "One or more supersaturation values to plot. "
            "If omitted, all values are shown with a shared colorbar."
        ),
    )
    p.add_argument(
        "--folders",
        "-f",
        nargs="+",
        default=None,
        metavar="FOLDER",
        help="Limit to specific subfolder names (default: all subfolders found).",
    )
    p.add_argument(
        "--rates",
        action="store_true",
        help=(
            "Also open a separate figure showing Δμ vs growth rate "
            "(linear slope of param over the x-axis) for each subfolder."
        ),
    )
    p.add_argument(
        "--rates-combined",
        action="store_true",
        help=(
            "Open a single combined Δμ vs growth rate figure with all subfolders "
            "overlaid, each as a distinct line/marker. Can be used with or without --rates."
        ),
    )
    p.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help=(
            "Save figures to file. Evolution figure uses this path; "
            "rates figure (if --rates) appends '_rates' before the extension."
        ),
    )
    p.add_argument(
        "--exclude-near-zero",
        "-e",
        nargs="*",
        type=float,
        default=None,
        metavar="BOUND",
        help=(
            "In evolution plots, drop data points whose parameter value falls inside an "
            "excluded band (and switch to a broken y-axis). With no value the band is "
            "(-1, 1). One positive value V excludes (0, V) (positive side only); one "
            "negative value V excludes (V, 0). Two values LO HI exclude (LO, HI), "
            "e.g. '-0.5 0.5' or '0.25 0.5'."
        ),
    )
    p.add_argument(
        "--exclude-rate-near-zero",
        "-E",
        nargs="*",
        type=float,
        default=None,
        metavar="BOUND",
        help=(
            "In rate plots, drop growth-rate points whose slope falls inside an excluded "
            "band. Same value semantics as --exclude-near-zero (no value -> (-1, 1); one "
            "value -> one side; two values -> LO HI)."
        ),
    )
    p.add_argument(
        "--symlog-x",
        action="store_true",
        help="Apply symmetric-log scale to the x-axis of all panels.",
    )
    p.add_argument(
        "--symlog-y",
        action="store_true",
        help="Apply symmetric-log scale to the y-axis of all panels.",
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="DPI for saved figures (default: 150).",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------


def find_csv_files(root: Path, folder_filter=None):
    """Return sorted list of (subfolder_name, csv_path) pairs."""
    results = []
    for csv_path in sorted(root.glob("**/RESULTS/size_analysis.csv")):
        subfolder = csv_path.parts[len(root.parts)]
        if folder_filter and subfolder not in folder_filter:
            continue
        results.append((subfolder, csv_path))
    return results


def load_data(csv_path: Path, param: str, supersat_filter=None):
    """Load and optionally filter a single size_analysis.csv."""
    df = pd.read_csv(csv_path)
    if supersat_filter is not None:
        df = df[df["x_supersat"].isin(supersat_filter)]
    df = df[df[param].notna()]
    return df


def moving_average(series: pd.Series, window: int) -> pd.Series:
    """Centred rolling mean; falls back to available data near edges."""
    return series.rolling(window=window, center=True, min_periods=1).mean()


def linear_slope(x: np.ndarray, y: np.ndarray) -> float:
    """Return the slope of the best-fit line through (x, y)."""
    if len(x) < 2:
        return np.nan
    slope, _ = np.polyfit(x, y, 1)
    return slope


def resolve_exclude_band(values):
    """Turn the CLI value list for an exclude filter into a (lo, hi) band.

    None       -> None          (filter disabled)
    []         -> (-1.0, 1.0)   (flag given without a value: default band)
    [V] V > 0  -> (0.0, V)      (exclude positive side only)
    [V] V < 0  -> (V, 0.0)      (exclude negative side only)
    [A, B]     -> (min, max)
    """
    if values is None:
        return None
    if len(values) == 0:
        return (-1.0, 1.0)
    if len(values) == 1:
        v = values[0]
        return (0.0, v) if v >= 0 else (v, 0.0)
    if len(values) == 2:
        return (min(values), max(values))
    sys.exit("--exclude-near-zero / --exclude-rate-near-zero take at most two values.")


# ---------------------------------------------------------------------------
# Colour helpers
# ---------------------------------------------------------------------------


def build_colormap(all_supersats):
    """Return (norm, cmap, sorted_supersats) consistent across all panels."""
    sorted_ss = sorted(all_supersats)
    vmin, vmax = sorted_ss[0], sorted_ss[-1]
    if vmin == vmax:
        vmin -= 0.5
        vmax += 0.5
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.coolwarm
    return norm, cmap, sorted_ss


def add_colorbar(fig, axes, cmap, norm):
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.02, pad=0.04)
    cbar.set_label(r"$\Delta\mu$", fontsize=12)
    return cbar


def add_colorbar_manual(fig, cmap, norm, rect=(0.90, 0.08, 0.02, 0.84)):
    """Colorbar at a fixed figure position — use when constrained_layout is off."""
    cax = fig.add_axes(rect)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="vertical")
    cbar.set_label(r"$\Delta\mu$", fontsize=12)
    return cbar


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------


def apply_symlog(ax, symlog_x: bool, symlog_y: bool):
    """Apply symmetric-log scale to x and/or y axes."""
    if symlog_x:
        ax.set_xscale("symlog")
    if symlog_y:
        ax.set_yscale("symlog")


def plot_evolution_panel(
    ax,
    df,
    param,
    xaxis,
    window,
    norm,
    cmap,
    sorted_supersats,
    symlog_x=False,
    symlog_y=False,
    exclude_near_zero=False,
):
    """Scatter raw data. If window is set, also draw a moving-average line per Δμ."""
    for ss in sorted_supersats:
        subset = df[df["x_supersat"] == ss].sort_values(xaxis)
        if exclude_near_zero:
            subset = subset[~subset[param].between(-1, 1, inclusive="both")]
        if subset.empty:
            continue
        color = cmap(norm(ss))
        x = subset[xaxis].values
        y = subset[param].values
        ax.scatter(x, y, color=color, s=6, alpha=0.35, linewidths=0)
        if window is not None:
            ax.plot(x, moving_average(pd.Series(y), window), color=color, linewidth=1.5)

    ax.set_xlabel(XAXIS_LABELS[xaxis], fontsize=9)
    ax.set_ylabel(PARAM_LABELS.get(param, param), fontsize=9)
    ax.tick_params(labelsize=8)
    apply_symlog(ax, symlog_x, symlog_y)


def plot_broken_evolution_panel(
    ax_top,
    ax_bot,
    df,
    param,
    xaxis,
    window,
    norm,
    cmap,
    sorted_supersats,
    symlog_x=False,
    symlog_y=False,
    exclude_band=(-1.0, 1.0),
):
    """
    Broken y-axis evolution panel.
    ax_top shows values > exclude_hi; ax_bot shows values < exclude_lo.
    The (exclude_lo, exclude_hi) band is excluded from both axes.
    """
    exclude_lo, exclude_hi = exclude_band
    top_ys, bot_ys = [], []

    for ss in sorted_supersats:
        subset = df[df["x_supersat"] == ss].sort_values(xaxis)
        subset = subset[~subset[param].between(exclude_lo, exclude_hi, inclusive="both")]
        if subset.empty:
            continue
        color = cmap(norm(ss))
        x = subset[xaxis].values
        y = subset[param].values

        top_mask = y > exclude_hi
        bot_mask = y < exclude_lo

        for ax, mask, ys_list in ((ax_top, top_mask, top_ys), (ax_bot, bot_mask, bot_ys)):
            if mask.any():
                ys_list.extend(y[mask].tolist())
                ax.scatter(x[mask], y[mask], color=color, s=6, alpha=0.35, linewidths=0)
                if window is not None:
                    ma = moving_average(pd.Series(y[mask]), window)
                    ax.plot(x[mask], ma, color=color, linewidth=1.5)

    # Set y-limits to the data range with 5 % padding, enforce the break gap
    def _set_ylim(ax, ys, side):
        if not ys:
            ax.set_axis_off()
            return
        lo, hi = min(ys), max(ys)
        pad = max((hi - lo) * 0.05, 0.5)
        if side == "top":
            ax.set_ylim(max(exclude_hi, lo - pad), hi + pad)
        else:
            ax.set_ylim(lo - pad, min(exclude_lo, hi + pad))

    _set_ylim(ax_top, top_ys, "top")
    _set_ylim(ax_bot, bot_ys, "bot")

    # Hide the inner spines and x-tick labels to create a clean break
    ax_top.spines["bottom"].set_visible(False)
    ax_bot.spines["top"].set_visible(False)
    ax_top.tick_params(labelbottom=False, bottom=False)

    add_break_marks(ax_top, ax_bot)

    ylabel = PARAM_LABELS.get(param, param)
    ax_top.set_ylabel(ylabel, fontsize=9)
    ax_bot.set_xlabel(XAXIS_LABELS[xaxis], fontsize=9)
    for ax in (ax_top, ax_bot):
        ax.tick_params(labelsize=8)
        if symlog_x:
            ax.set_xscale("symlog")
        if symlog_y:
            ax.set_yscale("symlog")


def plot_rate_panel(ax, df, param, xaxis, sorted_supersats, exclude_band=None):
    """
    Plot Δμ (x-axis) vs growth rate (y-axis) for one subfolder.
    Growth rate = linear slope of param vs xaxis for each supersat track.
    Points are coloured by their Δμ value using the same coolwarm scale.
    If exclude_band is set, growth rates inside (lo, hi) are dropped.
    """
    ss_vals, rates = [], []
    for ss in sorted_supersats:
        subset = df[df["x_supersat"] == ss].sort_values(xaxis)
        if len(subset) < 2:
            continue
        slope = linear_slope(subset[xaxis].values, subset[param].values)
        if exclude_band is not None and exclude_band[0] <= slope <= exclude_band[1]:
            continue
        ss_vals.append(ss)
        rates.append(slope)

    if not ss_vals:
        ax.text(
            0.5,
            0.5,
            "insufficient data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=8,
            color="grey",
        )
        return

    ss_arr = np.array(ss_vals)
    rate_arr = np.array(rates)

    # Colour each point by its own Δμ value
    vmin, vmax = ss_arr.min(), ss_arr.max()
    if vmin == vmax:
        vmin -= 0.5
        vmax += 0.5
    point_norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    colors = cm.coolwarm(point_norm(ss_arr))

    ax.scatter(ss_arr, rate_arr, c=colors, s=40, zorder=3, edgecolors="none")
    ax.plot(ss_arr, rate_arr, color="grey", linewidth=0.8, alpha=0.5, zorder=2)
    ax.axhline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.5)

    rate_label = f"d({PARAM_LABELS.get(param, param)}) / d({XAXIS_LABELS[xaxis]})"
    ax.set_xlabel(r"$\Delta\mu$", fontsize=10)
    ax.set_ylabel(rate_label, fontsize=9)
    ax.tick_params(labelsize=8)


def plot_combined_rates(ax, all_data, param, xaxis, sorted_supersats, symlog_y=False, exclude_band=None):
    """
    All subfolders' Δμ vs growth-rate curves on a single axes.
    Each subfolder gets a distinct colour and marker; Δμ is on the x-axis.
    If exclude_band is set, growth rates inside (lo, hi) are dropped.
    """
    # Cycle through tab10 colours and a set of markers
    colors = plt.get_cmap("tab10").colors
    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "+"]

    for i, (subfolder, df) in enumerate(all_data.items()):
        ss_vals, rates = [], []
        for ss in sorted_supersats:
            subset = df[df["x_supersat"] == ss].sort_values(xaxis)
            if len(subset) < 2:
                continue
            slope = linear_slope(subset[xaxis].values, subset[param].values)
            if exclude_band is not None and exclude_band[0] <= slope <= exclude_band[1]:
                continue
            ss_vals.append(ss)
            rates.append(slope)

        if not ss_vals:
            continue

        ss_arr = np.array(ss_vals)
        rate_arr = np.array(rates)
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]
        label = subfolder.replace("_", " ")

        # --- determine shift ---
        if len(rate_arr) > 1:
            # Option 2: robust default
            idx = np.argmin(np.abs(rate_arr))
            x0 = ss_arr[idx]
        else:
            x0 = 0.0

        # --- compute sigmoid-based shift ---
        # if len(ss_arr) >= 4:  # need enough points to fit
        #     x0 = fit_sigmoid_x0(ss_arr, rate_arr)
        # else:
        #     x0 = ss_arr[np.argmin(np.abs(rate_arr))]

        ss_shifted = ss_arr - x0

        ax.plot(ss_shifted, rate_arr, color=color, linewidth=1.2, alpha=0.7, zorder=2)
        ax.scatter(
            ss_shifted,
            rate_arr,
            color=color,
            marker=marker,
            s=45,
            zorder=3,
            edgecolors="none",
            label=label,
        )

        # smooth curve for visualisation
        # x_fit = np.linspace(ss_arr.min(), ss_arr.max(), 200)

        # try:
        #     popt, _ = curve_fit(sigmoid, ss_arr, rate_arr, maxfev=10000)
        #     y_fit = sigmoid(x_fit, *popt)

        #     ax.plot(x_fit - popt[2], y_fit, linestyle="--", color=color, alpha=0.6)
        # except ValueError:
        #     pass

    ax.axhline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.5)
    rate_label = f"d({PARAM_LABELS.get(param, param)}) / d({XAXIS_LABELS[xaxis]})"
    ax.set_xlabel(r"$\Delta\mu$", fontsize=11)
    ax.set_ylabel(rate_label, fontsize=10)
    ax.tick_params(labelsize=9)
    ax.legend(fontsize=9, framealpha=0.7)
    # if symlog_y:
    #     ax.set_yscale("symlog")


def make_figure(n_panels, panel_w=5, panel_h=4):
    ncols = min(2, n_panels)
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(panel_w * ncols, panel_h * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    return fig, axes, nrows, ncols


def make_broken_figure(n_panels, panel_w=5, panel_h=5):
    """Figure where every panel has a broken y-axis: top (y > 1) over bottom (y < -1)."""
    ncols = min(3, n_panels)
    nrows = (n_panels + ncols - 1) // ncols
    # Avoid constrained_layout — it collapses nested subgridspec axes.
    # Use explicit margins and spacing instead.
    fig = plt.figure(figsize=(panel_w * ncols, panel_h * nrows))
    outer_gs = fig.add_gridspec(
        nrows,
        ncols,
        hspace=0.45,
        wspace=0.35,
        left=0.08,
        right=0.88,
        top=0.92,
        bottom=0.06,
    )
    axes_top = np.empty((nrows, ncols), dtype=object)
    axes_bot = np.empty((nrows, ncols), dtype=object)
    for r in range(nrows):
        for c in range(ncols):
            inner = outer_gs[r, c].subgridspec(2, 1, hspace=0.05, height_ratios=[2, 1])
            axes_top[r, c] = fig.add_subplot(inner[0])
            axes_bot[r, c] = fig.add_subplot(inner[1], sharex=axes_top[r, c])
    return fig, axes_top, axes_bot, nrows, ncols


def add_break_marks(ax_top, ax_bot, d=0.012):
    """Draw diagonal slash marks between ax_top (bottom edge) and ax_bot (top edge)."""
    kw = dict(color="k", clip_on=False, linewidth=1, transform=ax_top.transAxes)
    ax_top.plot((-d, +d), (-d, +d), **kw)
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kw)
    kw["transform"] = ax_bot.transAxes
    ax_bot.plot((-d, +d), (1 - d, 1 + d), **kw)
    ax_bot.plot((1 - d, 1 + d), (1 - d, 1 + d), **kw)


def hide_unused(axes, n_used, nrows, ncols):
    for idx in range(n_used, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)


def hide_unused_broken(axes_top, axes_bot, n_used, nrows, ncols):
    for idx in range(n_used, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes_top[row][col].set_axis_off()
        axes_bot[row][col].set_axis_off()


def save_or_show(fig, output: Path, suffix: str, dpi: int):
    if output:
        stem = output.stem + suffix
        path = output.with_name(stem + output.suffix)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        print(f"Figure saved to {path}")
    else:
        fig.show()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    args = parse_args()

    if args.window is not None and args.window < 1:
        sys.exit("--window must be >= 1")

    csv_files = find_csv_files(args.root, folder_filter=args.folders)
    if not csv_files:
        sys.exit(f"No size_analysis.csv files found under {args.root}")

    # Load all data; collect full supersat universe
    all_supersats = set()
    raw_data = {}
    for subfolder, csv_path in csv_files:
        df = load_data(csv_path, args.param, supersat_filter=args.supersats)
        raw_data[subfolder] = df
        all_supersats.update(df["x_supersat"].unique())

    if not all_supersats:
        sys.exit("No data remaining after filtering. Check --supersats values.")

    norm, cmap, sorted_supersats = build_colormap(all_supersats)
    n = len(csv_files)

    evo_band = resolve_exclude_band(args.exclude_near_zero)
    rate_band = resolve_exclude_band(args.exclude_rate_near_zero)

    # ------------------------------------------------------------------ #
    # Figure 1 – time evolution                                           #
    # ------------------------------------------------------------------ #
    if evo_band is not None:
        fig_evo, axes_top, axes_bot, nrows_evo, ncols_evo = make_broken_figure(n)

        for idx, (subfolder, _) in enumerate(csv_files):
            row, col = divmod(idx, ncols_evo)
            plot_broken_evolution_panel(
                axes_top[row][col],
                axes_bot[row][col],
                raw_data[subfolder],
                args.param,
                args.xaxis,
                args.window,
                norm,
                cmap,
                sorted_supersats,
                symlog_x=args.symlog_x,
                symlog_y=args.symlog_y,
                exclude_band=evo_band,
            )
            axes_top[row][col].set_title(
                subfolder.replace("_", " "), fontsize=10, fontweight="bold"
            )

        hide_unused_broken(axes_top, axes_bot, n, nrows_evo, ncols_evo)
        add_colorbar_manual(fig_evo, cmap, norm)
        axes_evo = axes_top  # for suptitle reference only
    else:
        fig_evo, axes_evo, nrows_evo, ncols_evo = make_figure(n)

        for idx, (subfolder, _) in enumerate(csv_files):
            row, col = divmod(idx, ncols_evo)
            ax = axes_evo[row][col]
            plot_evolution_panel(
                ax,
                raw_data[subfolder],
                args.param,
                args.xaxis,
                args.window,
                norm,
                cmap,
                sorted_supersats,
                symlog_x=args.symlog_x,
                symlog_y=args.symlog_y,
                exclude_near_zero=False,
            )
            ax.set_title(subfolder.replace("_", " "), fontsize=10, fontweight="bold")

        hide_unused(axes_evo, n, nrows_evo, ncols_evo)
        add_colorbar(fig_evo, axes_evo, cmap, norm)

    title_parts = [f"{PARAM_LABELS.get(args.param, args.param)} vs {XAXIS_LABELS[args.xaxis]}"]
    if args.window is not None:
        title_parts.append(f"smoothing window = {args.window}")
    if args.supersats:
        title_parts.append(r"$\Delta\mu$ = " + ", ".join(str(s) for s in sorted(args.supersats)))
    fig_evo.suptitle("   |   ".join(title_parts), fontsize=12)

    # ------------------------------------------------------------------ #
    # Figure 2 – growth rates (optional)                                  #
    # ------------------------------------------------------------------ #
    if args.rates:
        fig_rate, axes_rate, nrows_rate, ncols_rate = make_figure(n)

        for idx, (subfolder, _) in enumerate(csv_files):
            row, col = divmod(idx, ncols_rate)
            ax = axes_rate[row][col]
            plot_rate_panel(
                ax,
                raw_data[subfolder],
                args.param,
                args.xaxis,
                sorted_supersats,
                exclude_band=rate_band,
            )
            ax.set_title(subfolder.replace("_", " "), fontsize=10, fontweight="bold")

        hide_unused(axes_rate, n, nrows_rate, ncols_rate)

        rate_title = (
            f"Growth Rate: d({PARAM_LABELS.get(args.param, args.param)}) "
            f"/ d({XAXIS_LABELS[args.xaxis]})   vs   $\\Delta\\mu$"
        )
        fig_rate.suptitle(rate_title, fontsize=12)

    # ------------------------------------------------------------------ #
    # Figure 3 – combined rates (optional)                               #
    # ------------------------------------------------------------------ #
    if args.rates_combined:
        fig_combined, ax_combined = plt.subplots(figsize=(7, 5), constrained_layout=True)
        plot_combined_rates(
            ax_combined,
            raw_data,
            args.param,
            args.xaxis,
            sorted_supersats,
            symlog_y=args.symlog_y,
            exclude_band=rate_band,
        )
        rate_title = (
            f"Growth Rate: d({PARAM_LABELS.get(args.param, args.param)}) "
            f"/ d({XAXIS_LABELS[args.xaxis]})   vs   $\\Delta\\mu$"
        )
        fig_combined.suptitle(rate_title, fontsize=12)

    # ------------------------------------------------------------------ #
    # Save or show                                                        #
    # ------------------------------------------------------------------ #
    if args.output:
        save_or_show(fig_evo, args.output, "", args.dpi)
        if args.rates:
            save_or_show(fig_rate, args.output, "_rates", args.dpi)
        if args.rates_combined:
            save_or_show(fig_combined, args.output, "_rates_combined", args.dpi)
    else:
        plt.show()


if __name__ == "__main__":
    main()
