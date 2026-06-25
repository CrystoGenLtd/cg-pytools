from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, Literal
import hashlib

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

from cgpytools.log import setup_logging

LOG = setup_logging(name="CG-PLOT")


# Configure Plot Themes
@dataclass
class PlotTheme:
    """Configuration for plot styling"""

    # Color schemes
    primary_colors: Dict[str, str] = None
    categorical_palette: list = None
    continuous_cmap: str = "viridis"
    diverging_cmap: str = "RdBu_r"

    # Typography
    font_family: str = "sans-serif"
    font_size_base: int = 11
    font_size_title: int = 14
    font_size_label: int = 12
    font_size_tick: int = 10
    font_weight_title: str = "bold"

    # Figure settings
    figure_dpi: int = 300
    figure_facecolor: str = "white"
    savefig_bbox: str = "tight"
    savefig_pad: float = 0.1

    # Grid and spines
    grid_alpha: float = 0.3
    grid_linewidth: float = 0.8
    spine_linewidth: float = 1.2

    # Markers and lines
    marker = "o"
    marker_alt = "s"
    marker_size: int = 10
    line_width: float = 2.0
    alpha_scatter: float = 0.7
    alpha_fill: float = 0.3

    def __post_init__(self):
        if self.primary_colors is None:
            self.primary_colors = {
                "blue": "#2E86AB",
                "red": "#F24236",
                "green": "#2E8B57",
                "orange": "#F18F01",
                "purple": "#8E44AD",
                "gray": "#7F8C8D",
                "teal": "#16A085",
                "pink": "#E91E63",
            }

        if self.categorical_palette is None:
            self.categorical_palette = [
                "#2E86AB",
                "#F24236",
                "#2E8B57",
                "#F18F01",
                "#8E44AD",
                "#16A085",
                "#E91E63",
                "#34495E",
            ]


class GlobalPlotStyler:
    """Global styling manager for all plots"""

    def __init__(
        self,
        theme: PlotTheme = None,
        style: Literal["modern", "classic", "minimal", "dark", "publication"] = "modern",
    ):
        """
        Initialize global plot styler

        Parameters
        ----------
        theme : PlotTheme, optional
            Custom theme configuration
        style : str, default="modern"
            Predefined style: "modern", "classic", "minimal", "dark", "publication"
        """
        self.theme = theme if theme else PlotTheme()
        self.style = style
        self._original_rcparams = None
        self._setup_style()

    def _setup_style(self):
        """Set up the global matplotlib style"""
        # Store original settings
        self._original_rcparams = mpl.rcParams.copy()

        # Apply base styling
        plt.style.use("default")  # Reset to clean slate

        # Apply predefined style
        style_configs = {
            "modern": self._modern_style,
            "classic": self._classic_style,
            "minimal": self._minimal_style,
            "dark": self._dark_style,
            "publication": self._publication_style,
        }

        if self.style in style_configs:
            style_configs[self.style]()
        else:
            self._modern_style()  # Default fallback

    def _modern_style(self):
        """Modern, clean style with subtle colors"""
        plt.rcParams.update(
            {
                # Figure
                "figure.facecolor": self.theme.figure_facecolor,
                "figure.edgecolor": "none",
                "figure.dpi": 100,  # Display DPI
                "savefig.dpi": self.theme.figure_dpi,
                "savefig.bbox": self.theme.savefig_bbox,
                "savefig.pad_inches": self.theme.savefig_pad,
                "savefig.facecolor": self.theme.figure_facecolor,
                # Font
                "font.family": self.theme.font_family,
                "font.size": self.theme.font_size_base,
                "axes.titlesize": self.theme.font_size_title,
                "axes.labelsize": self.theme.font_size_label,
                "xtick.labelsize": self.theme.font_size_tick,
                "ytick.labelsize": self.theme.font_size_tick,
                "legend.fontsize": self.theme.font_size_tick,
                "axes.titleweight": self.theme.font_weight_title,
                # Colors
                "axes.prop_cycle": plt.cycler("color", self.theme.categorical_palette),
                "axes.facecolor": "white",
                "axes.edgecolor": "#CCCCCC",
                "axes.linewidth": self.theme.spine_linewidth,
                # Grid
                "axes.grid": False,
                "axes.grid.axis": "both",
                "grid.alpha": self.theme.grid_alpha,
                "grid.linewidth": self.theme.grid_linewidth,
                "grid.color": "#E5E5E5",
                # Spines
                "axes.spines.top": False,
                "axes.spines.right": False,
                "axes.spines.left": True,
                "axes.spines.bottom": True,
                # Ticks
                "xtick.color": "#666666",
                "ytick.color": "#666666",
                "xtick.direction": "out",
                "ytick.direction": "out",
                "xtick.major.size": 4,
                "ytick.major.size": 4,
                "xtick.minor.size": 2,
                "ytick.minor.size": 2,
                # Legend
                "legend.frameon": True,
                "legend.framealpha": 0.9,
                "legend.facecolor": "white",
                "legend.edgecolor": "#CCCCCC",
                "legend.shadow": False,
            }
        )

    def _classic_style(self):
        """Traditional academic style"""
        self._modern_style()  # Start with modern
        plt.rcParams.update(
            {
                "axes.spines.top": True,
                "axes.spines.right": True,
                "axes.grid": False,
                "font.family": "serif",
                "mathtext.fontset": "cm",
            }
        )

    def _minimal_style(self):
        """Clean, minimal style"""
        self._modern_style()
        plt.rcParams.update(
            {
                "axes.grid": False,
                "axes.edgecolor": "#999999",
                "axes.linewidth": 0.8,
                "xtick.major.size": 2,
                "ytick.major.size": 2,
            }
        )

    def _dark_style(self):
        """Dark theme for presentations"""
        plt.rcParams.update(
            {
                "figure.facecolor": "#2E2E2E",
                "axes.facecolor": "#2E2E2E",
                "axes.edgecolor": "#CCCCCC",
                "axes.labelcolor": "white",
                "text.color": "white",
                "xtick.color": "white",
                "ytick.color": "white",
                "grid.color": "#555555",
                "savefig.facecolor": "#2E2E2E",
            }
        )
        self._modern_style()  # Apply other modern settings

    def _publication_style(self):
        """High-quality publication style"""
        self._modern_style()
        plt.rcParams.update(
            {
                "font.family": "serif",
                "font.size": 10,
                "axes.titlesize": 12,
                "axes.labelsize": 11,
                "xtick.labelsize": 9,
                "ytick.labelsize": 9,
                "legend.fontsize": 9,
                "lines.linewidth": 1.5,
                "axes.linewidth": 1.0,
                "figure.dpi": 150,
                "savefig.dpi": 600,  # High resolution for publications
            }
        )

    def apply_zingg_style(self, ax, title: str = "", show_legend: bool = False):
        """Apply consistent styling to Zingg plots"""
        # Classification lines
        ax.axhline(
            y=2 / 3,
            color=self.theme.primary_colors["blue"],
            linestyle="--",
            alpha=0.8,
            linewidth=1.5,
            label="2/3 threshold",
        )
        ax.axvline(
            x=2 / 3,
            color=self.theme.primary_colors["blue"],
            linestyle="--",
            alpha=0.8,
            linewidth=1.5,
        )

        # Styling
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_xlabel("S:M (Flatness)", fontweight="medium")
        ax.set_ylabel("M:L (Elongation)", fontweight="medium")

        if title:
            ax.set_title(title, pad=15)

        # Add shape classification regions as subtle background
        self._add_zingg_regions(ax)

        if show_legend:
            ax.legend(loc="upper right", framealpha=0.9)

    def _add_zingg_regions(self, ax):
        """Add subtle background regions for Zingg classification"""
        from matplotlib.patches import Rectangle

        alpha = 0.05

        # Define regions
        regions = {
            "Compact": Rectangle((0, 2 / 3), 2 / 3, 1 / 3, alpha=alpha, color="green"),
            "Platy": Rectangle((2 / 3, 2 / 3), 1 / 3, 1 / 3, alpha=alpha, color="blue"),
            "Bladed": Rectangle((0, 0), 2 / 3, 2 / 3, alpha=alpha, color="orange"),
            "Prolate": Rectangle((2 / 3, 0), 1 / 3, 2 / 3, alpha=alpha, color="red"),
        }

        for name, patch in regions.items():
            ax.add_patch(patch)

    def apply_correlation_style(self, ax, title: str = ""):
        """Apply styling to correlation heatmaps"""
        ax.set_title(title, pad=20)

        # Rotate x labels for better readability
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
        plt.setp(ax.get_yticklabels(), rotation=0)

    def get_color_palette(self, n_colors: int = None, palette_type: str = "categorical"):
        """Get appropriate color palette"""
        if palette_type == "categorical":
            if n_colors is None or n_colors <= len(self.theme.categorical_palette):
                return (
                    self.theme.categorical_palette[:n_colors]
                    if n_colors
                    else self.theme.categorical_palette
                )
            else:
                # Generate more colors using seaborn
                return sns.color_palette("husl", n_colors)
        elif palette_type == "continuous":
            return self.theme.continuous_cmap
        elif palette_type == "diverging":
            return self.theme.diverging_cmap
        else:
            return self.theme.categorical_palette

    def get_color_for_group(self, group: str, palette_type: str = "categorical") -> str:
        """
        Return a consistent color for a given group label.
        """
        palette = self.get_color_palette(palette_type=palette_type)
        # Stable hash via hashlib
        digest = hashlib.md5(group.encode("utf-8")).hexdigest()
        idx = int(digest, 16) % len(palette)
        return palette[idx]

    def save_figure(self, fig, filepath: Path, **kwargs):
        """Save figure with consistent settings"""
        save_kwargs = {
            "dpi": self.theme.figure_dpi,
            "bbox_inches": self.theme.savefig_bbox,
            "pad_inches": self.theme.savefig_pad,
            "facecolor": fig.get_facecolor(),
            "edgecolor": "none",
        }
        save_kwargs.update(kwargs)

        fig.savefig(filepath, **save_kwargs)

    def create_figure(
        self, figsize: Tuple[int, int] = (10, 8), **kwargs
    ) -> Tuple[plt.Figure, plt.Axes]:
        """Create figure with consistent styling"""
        fig_kwargs = {
            "figsize": figsize,
            "facecolor": self.theme.figure_facecolor,
            "edgecolor": "none",
        }
        fig_kwargs.update(kwargs)

        return plt.subplots(**fig_kwargs)

    def reset_style(self):
        """Reset to original matplotlib settings"""
        if self._original_rcparams:
            mpl.rcParams.update(self._original_rcparams)

    def __enter__(self):
        """Context manager entry"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - reset styling"""
        self.reset_style()


# Convenience function for quick setup
def setup_global_style(style: str = "modern", theme: PlotTheme = None) -> GlobalPlotStyler:
    """
    Quick setup for global plot styling

    Parameters
    ----------
    style : str
        Style name: "modern", "classic", "minimal", "dark", "publication"
    theme : PlotTheme, optional
        Custom theme configuration

    Returns
    -------
    GlobalPlotStyler
        Configured styler instance

    Examples
    --------
    >>> styler = setup_global_style("modern")
    >>> # All subsequent plots will use the modern style

    >>> # For temporary styling
    >>> with setup_global_style("publication") as styler:
    ...     # Create publication-quality plots here
    ...     pass
    """
    return GlobalPlotStyler(theme=theme, style=style)
