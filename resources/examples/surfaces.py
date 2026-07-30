"""
Example usage of :class:`cgpytools.analysis.surfaces.CrystalWulff`.

Run it from the repository root::

    python resources/examples/surfaces.py

To run the same operations on your own files, use the command line interface of
the module itself::

    python -m cgpytools.analysis.surfaces <json> info
    python -m cgpytools.analysis.surfaces <json> mesh -o shapes --reduced
"""

from pathlib import Path

import numpy as np

from cgpytools.analysis.surfaces import CrystalWulff, coeffs_to_xyz

REPO_ROOT = Path(__file__).resolve().parents[2]
JSON_PATH = REPO_ROOT / "resources" / "outputs" / "urea_surface_energies.json"
SAVEDIR = REPO_ROOT / "resources" / "outputs" / "example_wulff"


def example_usage(json_path: str | Path = JSON_PATH, savedir: str | Path = SAVEDIR):
    """
    Example of how to use the CrystalWulff class.

    Args:
        json_path: Path to a surface energies JSON file
        savedir: Directory the example Wulff mesh is written to
    """

    # Load the crystal structure along with any surface cut data in the file
    calculator = CrystalWulff.from_json(json_path)
    info = calculator.get_crystal_info()
    print(f"Loaded crystal with {info['n_symmetry_operations']} symmetry operations")

    # --- Path 1: facets and energies come from the JSON file itself ---
    if calculator.hkl is not None:
        print(f"Found {calculator.n_surfaces} surface cuts in {Path(json_path).name}")

        # Surface energies derived from the stored interaction energy counts
        surface_energies = calculator.facet_energies(calculator.energies["solvated"])
        calculator.log_surfaces(msg="example_usage", show_offset=True, console=True)

        # Wulff shape from the stored (solvated) energies, exported as an STL
        mesh = calculator.to_mesh(name="wulff_solvated", savedir=savedir)
        print(f"Wulff mesh: {mesh.vertices.shape[0]} vertices, volume {mesh.volume:.4f}")

        # The same shape as spherical harmonic coefficients, without saving
        coeffs = calculator.to_mesh(
            energies=surface_energies, save=False, to_return="sph", l_max=20
        )
        print(f"Spherical harmonic coefficients: {np.asarray(coeffs).shape}")

        # Reconstruct the shape from the harmonics and write it out as an xyz cloud
        calculator.to_xyz(name="wulff_solvated", savedir=savedir, l_max=20)
        points = coeffs_to_xyz(coeffs, write=False)
        print(f"Reconstructed point cloud: {points.shape}")

        # Keep only the lowest energy cut of each unique Miller index
        reduced = calculator.reduced_copy(context="solvated")
        print(f"Reduced from {calculator.n_surfaces} to {reduced.n_surfaces} unique facets")

        # The reduced instance behaves like any other calculator
        reduced_mesh = reduced.to_mesh(name="wulff_reduced", savedir=savedir)
        print(f"Reduced Wulff mesh volume: {reduced_mesh.volume:.4f}")

    # --- Path 2: facets and energies supplied externally (e.g. from a size file) ---
    external_facets = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 1, 1]])
    external_energies = np.array([1.2, 1.3, 1.1, 1.5, 1.8])

    external_mesh = calculator.calculate_wulff_shape(
        external_facets, external_energies, expand_symmetry=True, reduce_facets=True
    )
    print(f"Wulff mesh from external energies: volume {external_mesh.volume:.4f}")

    return calculator


if __name__ == "__main__":
    example_usage()
