"""cgpytools: crystal-shape / coarse-grained morphology analysis library.

Reusable modules for analysing crystal morphology from computational
simulations: shape characterisation, Wulff/surface processing, energy-network
parsing, and plotting utilities.
"""

from cgpytools.analysis.shape_analysis import ShapeAnalyser
from cgpytools.io.crystal import CrystalShape
from cgpytools.io.net import CGNet

__all__ = ["CGNet", "CrystalShape", "ShapeAnalyser"]

__version__ = "0.1.0"
