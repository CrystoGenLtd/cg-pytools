import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Union
from dataclasses import dataclass

from chmpy.crystal.wulff import WulffConstruction

from cgtools.log import setup_logging

LOG = setup_logging(name="WULFF")


@dataclass
class CrystalWulffCalculator:
    """
    Simplified class for crystallographic calculations and Wulff shape generation.
    Reads crystallographic data from JSON and works with external facet/energy data.
    """

    direct_matrix: np.ndarray
    reciprocal_matrix: np.ndarray
    symmetry_operations: np.ndarray

    @classmethod
    def from_json(cls, filename: Union[str, Path]):
        """
        Load crystallographic information from JSON file.

        Args:
            filename: Path to JSON file containing crystal structure data

        Returns:
            CrystalWulffCalculator instance
        """
        with Path(filename).open(encoding="utf-8") as f:
            j = json.load(f)

        # Support both formats from original class
        if "surface_cuts" in j:
            surface = j["surface_cuts"]["surface_energies"]["crystal"]
        elif "crystal" in j:
            surface = j["crystal"]
        else:
            surface = j["surface_energies"]["crystal"]

        # Extract space group symmetry operations
        sg = surface["space group"]
        symops = np.empty((len(sg["symmetry_operations"]), 3, 3))
        for i, sym in enumerate(sg["symmetry_operations"]):
            rots = np.array(sym["seitz"])[:3, :3]
            symops[i, :, :] = rots

        # Extract lattice matrices
        direct = np.array(surface["unit cell"]["direct_matrix"]).T
        reciprocal = np.array(surface["unit cell"]["reciprocal_matrix"]).T

        return cls(direct_matrix=direct, reciprocal_matrix=reciprocal, symmetry_operations=symops)

    def hkl_to_cart(self, planes: np.ndarray) -> np.ndarray:
        """
        Convert Miller indices to Cartesian coordinates.

        Args:
            planes: Array of Miller indices with shape (n_planes, 3)

        Returns:
            Normalized plane normal vectors in Cartesian coordinates
        """
        plane_normals_cart = planes @ self.reciprocal_matrix
        magnitudes = np.linalg.norm(plane_normals_cart, axis=1)
        plane_normals_cart /= magnitudes[:, np.newaxis]
        return plane_normals_cart

    def expand_facets(self, facets: np.ndarray) -> np.ndarray:
        """
        Expand facets using crystallographic symmetry operations.

        Args:
            facets: Array of Miller indices with shape (n_facets, 3)

        Returns:
            Expanded array of symmetry-equivalent facets
        """
        nfacets = facets.shape[0]
        nsymop = self.symmetry_operations.shape[0]
        expanded_facets = np.empty((nfacets * nsymop, 3), dtype=int)

        for i, s in enumerate(self.symmetry_operations):
            expanded_facets[i * nfacets : (i + 1) * nfacets] = (facets @ s.T).astype(int)

        return expanded_facets

    def reduce_facets(
        self, facets: np.ndarray, energies: np.ndarray, keep_negative: bool = True
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Reduce facets to unique crystallographic directions.

        Args:
            facets: Array of Miller indices with shape (n_facets, 3)
            energies: Corresponding surface energies
            keep_negative: Whether to keep negative Miller indices as separate facets

        Returns:
            Tuple of (reduced_facets, reduced_energies)
        """
        if facets.shape[0] != energies.shape[0]:
            raise ValueError("Number of facets and energies must match!")

        unique = {}

        for i, hkl in enumerate(facets):
            # Reduce to lowest common denominator
            gcd = np.gcd.reduce(np.abs(hkl))
            if gcd == 0:
                continue
            reduced = tuple(hkl // gcd)

            # Optionally consider negative directions as equivalent
            if not keep_negative:
                # Use positive form by convention
                if (
                    reduced[0] < 0
                    or (reduced[0] == 0 and reduced[1] < 0)
                    or (reduced[0] == 0 and reduced[1] == 0 and reduced[2] < 0)
                ):
                    reduced = tuple(-x for x in reduced)

            # Keep the facet with lowest energy for each unique direction
            if reduced not in unique or energies[i] < energies[unique[reduced]]:
                unique[reduced] = i

        # Extract unique facets and their energies
        indices = list(unique.values())
        reduced_facets = facets[indices]
        reduced_energies = energies[indices]

        LOG.debug(f"Reduced from {len(facets)} to {len(reduced_facets)} facets")
        return reduced_facets, reduced_energies

    def expand_and_reduce_facets(
        self, facets: np.ndarray, energies: np.ndarray, keep_negative: bool = True
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        First expand facets using symmetry operations, then reduce to unique directions.

        Args:
            facets: Array of Miller indices with shape (n_facets, 3)
            energies: Corresponding surface energies
            keep_negative: Whether to keep negative Miller indices as separate facets

        Returns:
            Tuple of (expanded_and_reduced_facets, corresponding_energies)
        """

        nsymop = self.symmetry_operations.shape[0]
        tiled_energies = np.tile(energies, nsymop)
        # Expand using symmetry
        expanded_facets = self.expand_facets(facets)
        # Reduce to unique directions
        return self.reduce_facets(expanded_facets, tiled_energies, keep_negative)

    def calculate_wulff_shape(
        self,
        facets: np.ndarray,
        energies: np.ndarray,
        expand_symmetry: bool = True,
        reduce_facets: bool = True,
        keep_negative: bool = True,
    ):  # -> Trimesh:  # Uncomment when Trimesh is available
        """
        Calculate Wulff shape from facets and surface energies.

        Args:
            facets: Array of Miller indices with shape (n_facets, 3)
            energies: Corresponding surface energies
            expand_symmetry: Whether to expand facets using symmetry operations
            reduce_facets: Whether to reduce to unique crystallographic directions
            keep_negative: Whether to keep negative Miller indices as separate facets

        Returns:
            Trimesh object representing the Wulff shape
        """
        # Process facets according to options
        if expand_symmetry:
            final_facets, final_energies = self.expand_and_reduce_facets(
                facets, energies, keep_negative
            )
        elif reduce_facets:
            final_facets, final_energies = self.reduce_facets(facets, energies, keep_negative)
        else:
            final_facets, final_energies = facets, energies

        # Convert to Cartesian coordinates
        plane_normals_cart = self.hkl_to_cart(final_facets)

        LOG.debug(f"Final facets shape: {final_facets.shape}")
        LOG.debug(f"Surface energies: {np.array2string(final_energies, precision=3)}")

        # Create Wulff construction
        wulff = WulffConstruction(
            facet_normals=plane_normals_cart,
            facet_energies=final_energies,
            labels=final_facets,
        )
        return wulff.to_trimesh()

    def get_crystal_info(self) -> dict:
        """
        Get summary of crystallographic information.

        Returns:
            Dictionary containing crystal structure information
        """
        return {
            "direct_matrix": self.direct_matrix,
            "reciprocal_matrix": self.reciprocal_matrix,
            "n_symmetry_operations": self.symmetry_operations.shape[0],
            "space_group_operations": self.symmetry_operations,
        }


def parse_miller_indices(hkl_string: str) -> np.ndarray:
    """
    Parse Miller indices from string format to numpy array.

    Args:
        hkl_string: String like "-1  0  0" or "0 -1  1"

    Returns:
        numpy array of Miller indices [h, k, l]
    """
    # Remove extra whitespace and split
    indices = hkl_string.strip().split()
    return np.array([int(idx) for idx in indices], dtype=int)


def parse_size_file_headers(size_file_path: str) -> np.ndarray:
    """
    Parse the Miller indices from the size file headers.

    Args:
        size_file_path: Path to the size file

    Returns:
        Array of Miller indices with shape (n_facets, 3)
    """
    with open(size_file_path, "r") as f:
        header_line = f.readline().strip()

    # Split by comma and remove 'time' column
    hkl_strings = [col.strip() for col in header_line.split(",")[1:]]

    facets = []
    for hkl_str in hkl_strings:
        if hkl_str:  # Skip empty strings
            facets.append(parse_miller_indices(hkl_str))

    return np.array(facets, dtype=int)


def generate_wulff_shapes_from_size_file(
    calculator: CrystalWulffCalculator,
    size_file_path: str | Path,
    n_steps: int = 10,
    file_format: str = "stl",
    expand_symmetry: bool = True,
    reduce_facets: bool = True,
    save_data: bool = True,
):
    """
    Generate Wulff shapes from size file data at specified intervals.

    Args:
        calculator: Instance of CrystalWulffCalculator class with the crystallographic information
        size_file_path: Path to the size file
        output_dir: Directory to save Wulff shapes
        n_steps: Generate Wulff shape every n time steps
        file_format: Output format ('obj', 'ply', 'stl', etc.)
        expand_symmetry: Whether to expand facets using crystal symmetry
        reduce_facets: Whether to reduce to unique crystallographic directions
        save_data: Whether to save intermediate data as CSV/JSON
    """

    # Read size file
    size_df = pd.read_csv(size_file_path)

    # Parse Miller indices from headers
    facets = parse_size_file_headers(size_file_path)

    # Create output directory
    size_dir = Path(size_file_path).parent
    output_path = size_dir / "shapes"
    output_path.mkdir(parents=True, exist_ok=True)

    LOG.info(f"Found {len(facets)} facets: {facets}")
    LOG.info(f"Processing {len(size_df)} time steps, saving every {n_steps} steps")

    # Store results for analysis
    results_data = []

    # Process every n-th time step
    for i in range(0, len(size_df), n_steps):
        timestep = i
        row = size_df.iloc[i]

        # Extract energies (skip the 'time' column)
        energies = row.iloc[1:].dropna().values.astype(float)

        # Ensure we have matching facets and energies
        if len(energies) != len(facets):
            LOG.warning(
                f"Mismatch at timestep {timestep}: {len(facets)} facets, {len(energies)} energies"
            )
            min_length = min(len(facets), len(energies))
            facets_subset = facets[:min_length]
            energies_subset = energies[:min_length]
        else:
            facets_subset = facets
            energies_subset = energies

        try:
            # Calculate Wulff shape
            wulff_result = calculator.calculate_wulff_shape(
                facets_subset,
                energies_subset,
                expand_symmetry=expand_symmetry,
                reduce_facets=reduce_facets,
            )

            # Save Wulff shape (when Trimesh is available)
            output_filename = f"wulff_timestep_{timestep:06d}.{file_format}"
            output_filepath = output_path / output_filename
            if output_filepath.exists():
                LOG.warning(f"Wulff shape file already exists, skipping: {output_filepath}")
            else:
                wulff_result.export(str(output_filepath))
                LOG.debug(f"Saved Wulff shape: {output_filepath}")

            result_info = {
                "timestep": timestep,
                "time_value": row.iloc[0],
                "n_facets": len(facets_subset),
                "n_energies": len(energies_subset),
                "output_file": str(output_filepath.absolute()),
            }
            # Add individual energy columns
            for i, energy in enumerate(energies_subset, 1):
                result_info[f"Int_{i}"] = energy

            results_data.append(result_info)

            LOG.debug(f"Processed timestep {timestep}: {len(facets_subset)} facets")

        except Exception as e:
            LOG.error(f"Error processing timestep {timestep}: {e}")
            continue

    # Save processing results
    if save_data and results_data:
        results_df = pd.DataFrame(
            [
                {
                    "timestep": r["timestep"],
                    "time_value": r["time_value"],
                    "n_facets": r["n_facets"],
                    "n_energies": r["n_energies"],
                    "output_file": r["output_file"],
                }
                for r in results_data
            ]
        )
        results_csv_path = output_path / "wulff_generation_summary.csv"
        results_df.to_csv(results_csv_path, index=False)

        # Save detailed data as JSON
        detailed_json_path = output_path / "wulff_detailed_data.json"
        with open(detailed_json_path, "w") as f:
            json.dump(results_data, f, indent=2, cls=NumpyEncoder)

        LOG.info(f"Saved summary: {results_csv_path}")
        LOG.info(f"Saved detailed data: {detailed_json_path}")

    return results_data


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder for numpy arrays."""

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        return super().default(obj)


def process_multiple_size_files(
    size_files: List[Path],
    crystal_json_path: str,
    base_output_dir: str,
    n_steps: int = 10,
    **kwargs,
):
    """
    Process multiple size files and generate Wulff shapes for each.

    Args:
        size_files: List of paths to size files
        crystal_json_path: Path to crystal structure JSON
        base_output_dir: Base directory for outputs
        n_steps: Generate Wulff shape every n time steps
        **kwargs: Additional arguments for generate_wulff_shapes_from_size_file
    """

    base_path = Path(base_output_dir)
    all_results = {}

    # Load crystal structure
    calculator = CrystalWulffCalculator.from_json(crystal_json_path)

    for size_file in size_files:
        if not isinstance(size_file, Path):
            size_file = Path(size_file)

        LOG.info(f"Processing size file: {size_file}")

        try:
            results = generate_wulff_shapes_from_size_file(
                calculator, size_file, n_steps=n_steps, **kwargs
            )
            all_results[str(size_file)] = results

            # Save all_results as CSV in base_path
            if all_results:
                combined_results = []
                for size_file_path, file_results in all_results.items():
                    size_file_name = Path(size_file_path).stem.replace("_size", "")
                    for result in file_results:
                        result_row = {
                            "size_file": size_file_name,
                            **result,  # Unpack all result data
                        }
                        combined_results.append(result_row)

                if combined_results:
                    combined_df = pd.DataFrame(combined_results)
                    csv_output_path = base_path / "wulff_processing_results.csv"
                    combined_df.to_csv(csv_output_path, index=False)
                    LOG.info(f"Saved combined Wulff processing results: {csv_output_path}")
                else:
                    LOG.warning("No valid results to save to CSV")

            else:
                LOG.warning("No results generated to save")

        except Exception as e:
            LOG.error(f"Failed to process {size_file}: {e}")
            continue

    return all_results


# Example usage:
def example_usage():
    """Example of how to use the CrystalWulffCalculator class."""

    # Load crystal structure
    calculator = CrystalWulffCalculator.from_json("crystal_data.json")

    # Define external facets and energies
    external_facets = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 1, 1]])

    external_energies = np.array([1.2, 1.3, 1.1, 1.5, 1.8])

    # Calculate Wulff shape
    calculator.calculate_wulff_shape(
        external_facets, external_energies, expand_symmetry=True, reduce_facets=True
    )

    # Get crystal info
    calculator.get_crystal_info()
    print("Crystal structure loaded successfully")

    return None


if __name__ == "__main__":
    example_usage()
