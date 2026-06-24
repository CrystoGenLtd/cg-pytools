"""cgtools: crystal-shape / coarse-grained morphology analysis library.

Reusable modules for analysing crystal morphology from computational
simulations: shape characterisation, Wulff/surface processing, energy-network
parsing, and plotting utilities.
"""

from cgtools.crystal_io import CrystalShape
from cgtools.shape_analysis import ShapeAnalyser
from cgtools.cg_net import CGNet

__all__ = ["CrystalShape", "ShapeAnalyser", "CGNet"]

__version__ = "0.1.0"
