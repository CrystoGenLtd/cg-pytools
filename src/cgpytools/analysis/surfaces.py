import argparse
import json
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, NamedTuple, NoReturn

import numpy as np
import pandas as pd
from chmpy.crystal.wulff import WulffConstruction
from chmpy.shape.reconstruct import reconstruct
from trimesh import Trimesh

from cgpytools.io.log import setup_logging

LOG = setup_logging(name="WULFF")

# eV/A^2 -> mJ/m^2 conversion used throughout the surface energy calculations
ENERGY_CONVERSION_FACTOR = 0.16604390679


def log_and_raise_error(msg: str, logger=LOG, exception=ValueError) -> NoReturn:
    """Log an error message then raise it as the given exception type."""
    logger.error(msg)
    raise exception(msg)


def normalise_verts(verts: np.ndarray, center: bool = True) -> np.ndarray:
    """
    Centre the vertices on their centroid and scale them onto the unit sphere.

    Args:
        verts: Cartesian coordinates with shape (n_points, 3)
        center: Whether to subtract the centroid before scaling

    Returns:
        The normalised vertices
    """
    verts = np.asarray(verts, dtype=float)
    if center:
        verts = verts - np.mean(verts, axis=0)
    norm = np.linalg.norm(verts, axis=1).max()
    return verts / norm


def coeffs_to_xyz(
    coeffs: np.ndarray,
    name: str = "wulff",
    index: int = 0,
    path: str | Path = ".",
    comment: str | None = None,
    suffix: str = "txt",
    normalise: bool = True,
    write: bool = True,
) -> np.ndarray:
    """
    Reconstruct a point cloud from spherical harmonic coefficients and write it out
    in xyz format.

    The file holds the number of points on the first line, a comment (the file path
    when none is given) on the second, and one ``x y z`` triple per line after that.

    Args:
        coeffs: Spherical harmonic coefficients
        name: Base name of the output file, written as ``{name}_{index}.{suffix}``
        index: Index appended to the filename, useful for trajectories
        path: Directory the file is written to
        comment: Comment line; defaults to the file path
        suffix: File extension for the output file
        normalise: Whether to centre and scale the reconstructed points
        write: Whether to write the file at all

    Returns:
        The reconstructed points with shape (n_points, 3)
    """
    LOG.debug("Reconstructing %s coeffs to xyz!", np.asarray(coeffs).shape)
    points = np.asarray(next(iter(reconstruct(coefficients=coeffs))), dtype=float)
    if normalise:
        points = normalise_verts(points)

    if write:
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        filepath = path / f"{name}_{index}.{suffix}"
        comment_line = comment if comment else str(filepath)
        with filepath.open("w", encoding="utf-8") as xyz_file:
            xyz_file.write(f"{len(points)}\n{comment_line}\n")
            for line in points:
                xyz_file.write(f"{line[0]}  {line[1]}  {line[2]}\n")

    return points


class FacetData(NamedTuple):
    """The facet arrays of a :class:`CrystalWulff`, known to be present."""

    hkl: np.ndarray
    offsets: np.ndarray
    area: np.ndarray
    facet_counts: np.ndarray


@dataclass
class CrystalWulff:
    """
    Simplified class for crystallographic calculations and Wulff shape generation.
    Reads crystallographic data from JSON and works with external facet/energy data.

    The facet data (``hkl``, ``offsets``, ``area``, ``facet_counts``, ``energies``)
    is optional: it is populated by :meth:`from_json` when the JSON file carries
    surface cut information, and is required by the interaction-energy based
    methods (:meth:`facet_energies`, :meth:`to_mesh`, ...).
    """

    direct_matrix: np.ndarray
    reciprocal_matrix: np.ndarray
    symmetry_operations: np.ndarray
    hkl: np.ndarray | None = None
    offsets: np.ndarray | None = None
    area: np.ndarray | None = None
    facet_counts: np.ndarray | None = None
    energies: dict[str, np.ndarray] = field(default_factory=dict)

    @classmethod
    def from_json(cls, filename: str | Path):
        """
        Load crystallographic information from JSON file.

        Args:
            filename: Path to JSON file containing crystal structure data

        Returns:
            CrystalWulff instance
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

        hkl, offsets, area, facet_counts, energies = cls._parse_surface_cuts(j)

        return cls(
            direct_matrix=direct,
            reciprocal_matrix=reciprocal,
            symmetry_operations=symops,
            hkl=hkl,
            offsets=offsets,
            area=area,
            facet_counts=facet_counts,
            energies=energies,
        )

    @staticmethod
    def _parse_surface_cuts(j: dict):
        """
        Parse the facet/interaction energy blocks of a surface energies JSON file.

        Returns a tuple of ``(hkl, offsets, area, facet_counts, energies)``, all of
        which are ``None``/empty when the file carries no surface cut information.
        """

        def to_arr_with_consistent_shape(list_of_dicts):
            try:
                list_of_values = [[x["total"] for x in l] for l in list_of_dicts]
            except KeyError:
                list_of_values = [[x["Total"] for x in l] for l in list_of_dicts]
            except TypeError:
                # Old format: already a list of lists of floats
                list_of_values = list_of_dicts

            arrs = [np.trim_zeros(x) for x in list_of_values]
            dims = len(arrs), max(len(arr) for arr in arrs)
            result = np.zeros(dims)
            for i, arr in enumerate(arrs):
                result[i, : len(arr)] = arr
            return result

        if "surface_cuts" in j:
            cuts = j["surface_cuts"]
            facets = cuts["surface_energies"]["facets"]
        elif "vacuum" in j and "surface_energies" in j:
            cuts = j
            facets = j["surface_energies"]["facets"]
        else:
            LOG.debug("No surface cut data found in JSON; facet data unavailable")
            return None, None, None, None, {}

        vacuum_energies = -to_arr_with_consistent_shape(cuts["vacuum"])
        solvated_energies = to_arr_with_consistent_shape(cuts["solvated"])

        num_facets = len(facets)
        hkl = np.empty((num_facets, 3), dtype=int)
        offsets = np.empty(num_facets)
        area = np.empty(num_facets)
        facet_counts = np.zeros((num_facets, *vacuum_energies.shape))

        for i, facet in enumerate(facets):
            hkl[i] = facet["hkl"]
            offsets[i] = facet["offset"]
            area[i] = facet["area"]
            for k, row in enumerate(facet["interaction_energy_counts"]):
                n = min(len(row), facet_counts.shape[2])
                facet_counts[i, k, :n] = row[:n]

        energies = {"vacuum": vacuum_energies, "solvated": solvated_energies}
        LOG.debug("Parsed %s surface cuts from JSON", num_facets)

        return hkl, offsets, area, facet_counts, energies

    @property
    def n_surfaces(self) -> int:
        return self._require_facet_data().hkl.shape[0]

    @property
    def n_molecules(self) -> int:
        self._require_facet_data()
        return self.energies["vacuum"].shape[0]

    @property
    def n_energies(self) -> int:
        self._require_facet_data()
        return self.energies["vacuum"].shape[1]

    def _require_facet_data(self) -> FacetData:
        """
        Return the facet arrays, raising if the instance carries no facet/interaction
        energy information.
        """
        if (
            self.hkl is None
            or self.offsets is None
            or self.area is None
            or self.facet_counts is None
            or not self.energies
        ):
            log_and_raise_error(
                "No facet data available: this instance was created without surface "
                "cut information (use from_json on a surface energies JSON file).",
                LOG,
                ValueError,
            )
        return FacetData(self.hkl, self.offsets, self.area, self.facet_counts)

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

    def expand_facets(self, facets: np.ndarray | None = None) -> np.ndarray:
        """
        Expand facets using crystallographic symmetry operations.

        Args:
            facets: Array of Miller indices with shape (n_facets, 3).
                Defaults to the facets loaded from JSON (``self.hkl``).

        Returns:
            Expanded array of symmetry-equivalent facets
        """
        if facets is None:
            facets = self._require_facet_data().hkl

        nfacets = facets.shape[0]
        nsymop = self.symmetry_operations.shape[0]
        expanded_facets = np.empty((nfacets * nsymop, 3), dtype=int)

        for i, s in enumerate(self.symmetry_operations):
            expanded_facets[i * nfacets : (i + 1) * nfacets] = (facets @ s.T).astype(int)

        return expanded_facets

    def reduce_facets(
        self, facets: np.ndarray, energies: np.ndarray, keep_negative: bool = True
    ) -> tuple[np.ndarray, np.ndarray]:
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
            # Use positive form by convention
            if not keep_negative and (
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

        LOG.debug("Reduced from %s to %s facets", len(facets), len(reduced_facets))
        return reduced_facets, reduced_energies

    def expand_and_reduce_facets(
        self, facets: np.ndarray, energies: np.ndarray, keep_negative: bool = True
    ) -> tuple[np.ndarray, np.ndarray]:
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

        LOG.debug("Final facets shape: %s", final_facets.shape)
        LOG.debug("Surface energies: %s", np.array2string(final_energies, precision=3))

        # Create Wulff construction
        wulff = WulffConstruction(
            facet_normals=plane_normals_cart,
            facet_energies=final_energies,
            labels=final_facets,
        )
        return wulff.to_trimesh()

    def facet_energies(self, interaction_energies: np.ndarray) -> np.ndarray:
        """
        Convert per-molecule interaction energy counts into surface energies.

        Args:
            interaction_energies: Interaction energies with shape
                (n_molecules, n_energies)

        Returns:
            Surface energy per facet, shape (n_surfaces,)
        """
        _, _, area, facet_counts = self._require_facet_data()
        try:
            etot = np.einsum("ijk,jk->i", facet_counts, interaction_energies)
        except ValueError as exc:
            log_and_raise_error(
                "Couldn't perform Einstein summation: "
                f"{exc}\n Facet counts: {facet_counts.shape} | "
                f"Interaction energies: {np.asarray(interaction_energies).shape}",
                LOG,
                ValueError,
            )
        return ENERGY_CONVERSION_FACTOR * 0.5 * etot / area

    def expand_symmetry_related_planes(self, energies: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Expand the loaded facets by symmetry, keeping the lowest energy for each
        unique crystallographic direction.

        Args:
            energies: Surface energies, one per loaded facet (shape (n_surfaces,))

        Returns:
            Tuple of (unique_facets, unique_energies)
        """
        energies = np.asarray(energies)
        if energies.shape != (self.n_surfaces,):
            LOG.error(
                "Number of facets found: %s | Number of energies received: %s",
                self.n_surfaces,
                energies.shape,
            )
            raise ValueError("Number of facets and number of energies received does not match!")

        expanded_facets = self.expand_facets()
        nsymop = self.symmetry_operations.shape[0]
        tiled_energies = np.tile(energies, nsymop)

        # identify the unique directions
        unique = {}
        try:
            for i, x in enumerate(expanded_facets):
                reduced = tuple(x / np.gcd.reduce(x).max())
                if reduced not in unique or tiled_energies[unique[reduced]] > tiled_energies[i]:
                    unique[reduced] = i
        except IndexError as ie:
            raise IndexError(f"with {tiled_energies.shape=} | {expanded_facets.shape=}") from ie

        unique_facets = []
        unique_energies = []
        for k, v in unique.items():
            unique_facets.append(k)
            unique_energies.append(tiled_energies[v])

        return np.array(unique_facets), np.array(unique_energies)

    def wulff_shape_sht(self, planes: np.ndarray, energies: np.ndarray, l_max: int = 20):
        """
        Spherical harmonic expansion of the Wulff shape for the given planes.

        Args:
            planes: Miller indices of the (already expanded) facets
            energies: Surface energy of each plane
            l_max: Maximum degree of the expansion

        Returns:
            Spherical harmonic coefficients
        """
        plane_normals_cart = self.hkl_to_cart(planes)
        LOG.debug("Surface Energies: %s", np.array2string(energies, precision=2))
        wulff = WulffConstruction(
            facet_normals=plane_normals_cart,
            facet_energies=energies,
            labels=planes,
        )
        return wulff.sht(l_max=l_max).coeffs

    def wulff_shape_mesh(self, planes: np.ndarray, energies: np.ndarray) -> Trimesh:
        """
        Wulff shape mesh for the given planes.

        Args:
            planes: Miller indices of the (already expanded) facets
            energies: Surface energy of each plane

        Returns:
            Trimesh of the Wulff shape
        """
        plane_normals_cart = self.hkl_to_cart(planes)
        LOG.debug("Surface Energies: %s", np.array2string(energies, precision=2))
        wulff = WulffConstruction(
            facet_normals=plane_normals_cart,
            facet_energies=energies,
            labels=planes,
        )
        return wulff.to_trimesh()

    def interaction_energies_to_planes(
        self, interaction_energies: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Convert interaction energies to symmetry expanded planes and energies."""
        energies = self.facet_energies(interaction_energies)
        return self.expand_symmetry_related_planes(energies)

    def surface_energies_to_coeffs(self, facet_energies: np.ndarray, l_max: int = 20):
        """Spherical harmonic coefficients of the Wulff shape from surface energies."""
        planes, energies = self.expand_symmetry_related_planes(facet_energies)
        return self.wulff_shape_sht(planes, energies, l_max=l_max)

    def interaction_energies_to_coeffs(self, interaction_energies: np.ndarray, l_max: int = 20):
        """Spherical harmonic coefficients of the Wulff shape from interaction energies."""
        planes, energies = self.interaction_energies_to_planes(interaction_energies)
        return self.wulff_shape_sht(planes, energies, l_max=l_max)

    def any_energies_to_surface_energies(
        self, energies: np.ndarray | list, surface: bool | None = None
    ) -> np.ndarray:
        """
        Ensure the provided energies are surface energies, converting interaction
        energies with :meth:`facet_energies` when required.

        Args:
            energies: Surface or interaction energies
            surface: Set explicitly when the two energy arrays have the same size
                and the type cannot be inferred

        Returns:
            Surface energies, one per facet
        """
        checked_energies, surface_detected = self.check_surface_vs_interaction_energies(
            energies, surface
        )

        if not surface_detected:
            return self.facet_energies(checked_energies)

        return checked_energies

    def check_surface_vs_interaction_energies(
        self, energies: np.ndarray | list, surface: bool | None = None
    ) -> tuple[np.ndarray, bool]:
        """
        Determine whether the given energies are surface or interaction energies,
        reshaping them to the expected shape.

        Returns:
            Tuple of (reshaped_energies, is_surface_energies)
        """
        self._require_facet_data()
        energies_arr = np.array(energies, dtype=float)
        interaction_energies = self.energies["vacuum"]
        surface_energies = self.facet_energies(interaction_energies)

        n_interaction_energies = interaction_energies.size
        n_surface_energies = surface_energies.size

        if n_interaction_energies == n_surface_energies:
            if surface is None:
                log_and_raise_error(
                    "Ambiguity detected: Please specify 'surface=True' or 'surface=False'.",
                    LOG,
                    ValueError,
                )

            correct_size = n_surface_energies if surface else n_interaction_energies
            if energies_arr.size != correct_size:
                log_and_raise_error(
                    "Mismatch in the size of provided "
                    f"{'surface' if surface else 'interaction'} energies.",
                    LOG,
                    ValueError,
                )

        else:
            if energies_arr.size == n_interaction_energies:
                surface = False
                LOG.debug("Energies Check (Passed): interaction energies")
            elif energies_arr.size == n_surface_energies:
                surface = True
                LOG.debug("Energies Check (Passed): surface energies")
            else:
                log_and_raise_error(
                    f"Provided energies {energies_arr.size} do not match known interaction "
                    f"({n_interaction_energies}) or surface ({n_surface_energies}) energy sizes.",
                    LOG,
                    ValueError,
                )

        # Reshape energies to match the correct shape before returning
        target_shape = surface_energies.shape if surface else interaction_energies.shape
        if energies_arr.shape != target_shape:
            LOG.info("Reshaping energies from %s to %s", energies_arr.shape, target_shape)
            energies_arr = np.reshape(energies_arr, target_shape)

        return energies_arr, surface

    def to_mesh(
        self,
        energies: np.ndarray | None = None,
        name: str = "wulff",
        savedir: Path | str = ".",
        suffix: str = ".stl",
        surface: bool | None = None,
        l_max: int = 20,
        solvation: Literal["vacuum", "solvated"] = "solvated",
        save: bool = True,
        *,
        to_return: Literal["mesh", "sph"] = "mesh",
    ) -> Trimesh | np.ndarray:
        """
        Generate the Wulff shape from the loaded facets and export it to a file.

        Args:
            energies: Surface or interaction energies; defaults to
                ``self.energies[solvation]``
            name: Base name for the output file
            savedir: Directory the output file is written to
            suffix: File suffix, determines the export format
            surface: Whether the given energies are surface energies (only needed
                when the type cannot be inferred from the array size)
            l_max: Maximum degree for the spherical harmonic expansion
            solvation: Which stored energies to use when ``energies`` is None
            save: Whether to export the mesh to disk
            to_return: Return the Trimesh ("mesh") or the SHT coefficients ("sph")

        Returns:
            The Wulff Trimesh, or the spherical harmonic coefficients
        """
        hkl = self._require_facet_data().hkl

        if energies is None:
            energies = self.energies[solvation]

        savedir = Path(savedir)

        _energies = self.any_energies_to_surface_energies(energies, surface)

        self.log_surfaces(hkl, _energies, "to_mesh")

        equivalent_planes, facet_energies = self.expand_symmetry_related_planes(_energies)

        plane_normals_cart = self.hkl_to_cart(equivalent_planes)

        wulff = WulffConstruction(
            facet_normals=plane_normals_cart,
            facet_energies=facet_energies,
            labels=equivalent_planes,
        )
        wulff_mesh = wulff.to_trimesh()

        if save:
            savedir.mkdir(parents=True, exist_ok=True)

            path = savedir / f"{name}{suffix}"
            wulff_mesh.export(path)
            LOG.debug("Saving %s file %s to: %s", suffix, name, savedir)

        if to_return == "mesh":
            return wulff_mesh
        elif to_return == "sph":
            return wulff.sht(l_max=l_max).coeffs
        else:
            raise ValueError(f"Specify 'to_return' as 'sph' or 'mesh'. Received {to_return}")

    def to_xyz(
        self,
        energies: np.ndarray | None = None,
        name: str = "wulff",
        comment: str | None = None,
        savedir: Path | str = ".",
        suffix: str = "txt",
        surface: bool | None = None,
        l_max: int = 20,
        solvation: Literal["vacuum", "solvated"] = "solvated",
        index: int = 0,
        normalise: bool = True,
        write: bool = True,
    ) -> np.ndarray:
        """
        Convert energies to spherical harmonic coefficients and save the reconstructed
        point cloud in xyz format.

        The energies are first classified as surface or interaction energies, expanded
        by symmetry, transformed into spherical harmonic coefficients and finally
        reconstructed into Cartesian points which are written to disk.

        Args:
            energies: Surface or interaction energies; defaults to
                ``self.energies[solvation]``
            name: Base name for the output file (written as ``{name}_{index}.{suffix}``)
            comment: Comment line written into the file; defaults to the file path
            savedir: Directory the output file is written to
            suffix: File extension for the output file
            surface: Whether the given energies are surface energies (only needed when
                the type cannot be inferred from the array size)
            l_max: Maximum degree of the spherical harmonic expansion
            solvation: Which stored energies to use when ``energies`` is None
            index: Index appended to the filename, useful for trajectories
            normalise: Whether to centre and scale the reconstructed points
            write: Whether to write the xyz file

        Returns:
            The spherical harmonic coefficients used to build the point cloud
        """
        self._require_facet_data()

        if energies is None:
            energies = self.energies[solvation]

        # Check and label energy types
        energies, surface_detected = self.check_surface_vs_interaction_energies(energies, surface)

        savedir = Path(savedir)

        # Convert energies to spherical harmonic coefficients
        if surface_detected:
            coeffs = self.surface_energies_to_coeffs(energies, l_max=l_max)
        else:
            coeffs = self.interaction_energies_to_coeffs(energies, l_max=l_max)

        coeffs_to_xyz(
            coeffs,
            name=name,
            index=index,
            path=savedir,
            comment=comment,
            suffix=suffix,
            normalise=normalise,
            write=write,
        )
        if write:
            LOG.debug("Saving .%s file %s to: %s", suffix, name, savedir)

        return coeffs

    def _check_keep(self, keep: list[tuple] | None = None) -> set[tuple]:
        """
        Validate the requested subset of facets, warning about any hkl that is not
        present on this instance.

        Args:
            keep: Miller indices to keep; all loaded facets are kept if None

        Returns:
            Set of hkl tuples to keep
        """
        _all = {tuple(hkl) for hkl in self._require_facet_data().hkl}
        if keep is None:
            return _all

        _keep = {tuple(hkl) for hkl in keep if tuple(hkl) in _all}
        if len(keep) != len(_keep):
            LOG.warning(
                "Only %s facets found! Found: %s Received (%s): %s",
                len(_keep),
                _keep,
                len(keep),
                keep,
            )
        return _keep

    def reduce_hkl(self, keep: list[tuple] | None = None) -> np.ndarray:
        """Return the unique Miller indices of the loaded facets."""
        kept = self._check_keep(keep)

        unique_hkl_dict = {}
        for hkl in self._require_facet_data().hkl:
            hkl_tuple = tuple(hkl)
            if hkl_tuple in kept:
                unique_hkl_dict[hkl_tuple] = []

        return np.array(list(unique_hkl_dict.keys()), dtype=np.int64)

    def get_min_energies(self, energies: np.ndarray, keep: list[tuple] | None = None) -> dict:
        """
        Compute the minimum energy cut for each unique Miller index.

        Args:
            energies: Surface energy of each loaded facet
            keep: Miller indices to consider; all facets are used if None

        Returns:
            Dictionary with the unique ``hkl``, the ``indices`` of the minimum energy
            cuts and the corresponding ``energies``
        """
        kept = self._check_keep(keep)

        min_energies_dict = OrderedDict()

        # Determine the minimum energy for each unique hkl
        for i, hkl in enumerate(self._require_facet_data().hkl):
            hkl_tuple = tuple(hkl)
            if hkl_tuple not in kept:
                continue
            if hkl_tuple not in min_energies_dict or energies[i] < min_energies_dict[hkl_tuple][1]:
                min_energies_dict[hkl_tuple] = (i, energies[i])

        sorted_hkl = np.array(list(min_energies_dict.keys()), dtype=np.int64)
        indices_and_energies = np.array(list(min_energies_dict.values()), dtype=np.float64)
        indices, min_energies = indices_and_energies[:, 0], indices_and_energies[:, 1]

        return {
            "hkl": sorted_hkl,
            "indices": np.int64(indices),
            "energies": min_energies,
        }

    def reduced_surface_energies(
        self,
        return_min_index: bool = False,
        context: Literal["vacuum", "solvated"] = "solvated",
        keep: list[tuple] | None = None,
    ) -> dict:
        """
        Minimum surface energy per unique Miller index, in both vacuum and solvated
        contexts.

        Args:
            return_min_index: Whether to include the indices of the minimum energy cuts
            context: Which energies determine the minimum energy cut
            keep: Miller indices to consider; all facets are used if None

        Returns:
            Dictionary with the ``vacuum`` and ``solvated`` energies, the selected
            ``hkl`` and, optionally, the ``min_indices``
        """
        energies_vacuum = self.facet_energies(self.energies["vacuum"])
        energies_solvated = self.facet_energies(self.energies["solvated"])

        min_energies_vacuum = self.get_min_energies(np.array(energies_vacuum), keep)
        min_energies_solvated = self.get_min_energies(np.array(energies_solvated), keep)

        if context == "solvated":
            selected = min_energies_solvated
        elif context == "vacuum":
            selected = min_energies_vacuum
        else:
            log_and_raise_error(
                f"Invalid context provided: {context}. Expected 'vacuum' or 'solvated'.",
                LOG,
                ValueError,
            )

        self.log_surfaces(
            selected["hkl"],
            selected["energies"],
            f"Unique surfaces ({context})",
            console=False,
            show_offset=False,
        )

        return {
            "vacuum": min_energies_vacuum["energies"],
            "solvated": min_energies_solvated["energies"],
            "min_indices": selected["indices"] if return_min_index else None,
            "hkl": selected["hkl"],
        }

    def filter_attributes(
        self,
        min_index: np.ndarray | None = None,
        keep: list[tuple] | None = None,
    ) -> OrderedDict:
        """
        Collect the offset, area and interaction energy counts of the selected cuts,
        keyed by their Miller index.

        Args:
            min_index: Indices of the cuts to keep
            keep: Miller indices to consider; all facets are used if None

        Returns:
            Ordered mapping of hkl tuple -> [offset, area, facet_counts]
        """
        kept = self._check_keep(keep)
        facet_data = self._require_facet_data()

        filtered_values_dict: OrderedDict[tuple, list | None] = OrderedDict(
            [(tuple(hkl), None) for hkl in facet_data.hkl if tuple(hkl) in kept]
        )

        for i, hkl in enumerate(facet_data.hkl):
            hkl_tuple = tuple(hkl)
            if min_index is not None and i in min_index:
                filtered_values_dict[hkl_tuple] = [
                    facet_data.offsets[i],
                    facet_data.area[i],
                    facet_data.facet_counts[i],
                ]

        return filtered_values_dict

    def reduced_copy(
        self,
        mode: Literal["min"] = "min",
        context: Literal["vacuum", "solvated"] = "solvated",
        keep: list[tuple] | None = None,
    ) -> "CrystalWulff":
        """
        Create a copy of this instance holding only one cut per unique Miller index.

        Args:
            mode: Selection rule; only "min" (lowest surface energy) is supported
            context: Which energies determine the minimum energy cut
            keep: Miller indices to keep; all facets are kept if None

        Returns:
            A new CrystalWulff with the reduced facet data
        """
        if mode != "min":
            raise NotImplementedError(
                "Currently only allowed to reduce based on min surface energies!"
            )

        min_surface_energies = self.reduced_surface_energies(
            return_min_index=True, context=context, keep=keep
        )
        reduced_hkl = min_surface_energies["hkl"]
        LOG.debug("Reduced HKL: %s", reduced_hkl)
        indices = min_surface_energies["min_indices"]

        surface_attributes = list(self.filter_attributes(indices, keep=keep).values())
        facet_data = self._require_facet_data()
        n = len(reduced_hkl)

        reduced_offsets = np.empty(n)
        reduced_area = np.empty(n)
        reduced_facet_counts = np.empty(shape=(n, self.n_molecules, self.n_energies))

        for i, x in enumerate(surface_attributes):
            if x is None:
                continue
            reduced_offsets[i] = x[0]
            reduced_area[i] = x[1]
            reduced_facet_counts[i] = np.array(x[2])

        LOG.debug(
            "%s -> %s facets (offsets %s -> %s, area %s -> %s, counts %s -> %s)",
            facet_data.hkl.shape,
            reduced_hkl.shape,
            facet_data.offsets.shape,
            reduced_offsets.shape,
            facet_data.area.shape,
            reduced_area.shape,
            facet_data.facet_counts.shape,
            reduced_facet_counts.shape,
        )

        return CrystalWulff(
            direct_matrix=self.direct_matrix,
            reciprocal_matrix=self.reciprocal_matrix,
            symmetry_operations=self.symmetry_operations,
            hkl=reduced_hkl,
            offsets=reduced_offsets,
            area=reduced_area,
            facet_counts=reduced_facet_counts,
            energies=self.energies,
        )

    def log_surfaces(
        self,
        hkls: np.ndarray | None = None,
        energies: np.ndarray | None = None,
        msg: str | None = None,
        solvation: Literal["vacuum", "solvated"] = "solvated",
        show_offset: bool = False,
        console: bool = False,
    ):
        """
        Log the surface planes and their energies, optionally printing to console.

        Args:
            hkls: Miller indices; uses ``self.hkl`` if None
            energies: Corresponding surface energies; computed from
                ``self.energies[solvation]`` if None
            msg: Additional message to prepend to the logs
            solvation: Context of the energies used when ``energies`` is None
            show_offset: Whether to include the facet offset in the log
            console: If True, prints the information as well as logging it
        """
        facet_data = self._require_facet_data()

        hkls = hkls if hkls is not None else facet_data.hkl
        energies = (
            energies if energies is not None else self.facet_energies(self.energies[solvation])
        )
        energies = np.asarray(energies).flatten()

        if msg is not None:
            log_msg = f"Surfaces from: {msg}"
            LOG.debug(log_msg)
            if console:
                print(log_msg)

        log_msg = f"[HKL] {hkls.shape} -> ENERGY {energies.shape}"
        LOG.debug(log_msg)
        if console:
            print(log_msg)

        order = np.argsort(energies)
        for idx in order:
            log_msg = f"{np.array2string(hkls[idx]):>12}  ->  {energies[idx]}   "
            log_msg += f"{facet_data.offsets[idx] if show_offset else ''}"

            LOG.debug(log_msg)
            if console:
                print(log_msg)

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


def parse_size_file_headers(size_file_path: str | Path) -> np.ndarray:
    """
    Parse the Miller indices from the size file headers.

    Args:
        size_file_path: Path to the size file

    Returns:
        Array of Miller indices with shape (n_facets, 3)
    """
    with open(size_file_path, "r", encoding="utf-8") as f:
        header_line = f.readline().strip()

    # Split by comma and remove 'time' column
    hkl_strings = [col.strip() for col in header_line.split(",")[1:]]

    facets = []
    for hkl_str in hkl_strings:
        if hkl_str:  # Skip empty strings
            facets.append(parse_miller_indices(hkl_str))

    return np.array(facets, dtype=int)


def generate_wulff_shapes_from_size_file(
    calculator: CrystalWulff,
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
        calculator: Instance of CrystalWulff class with the crystallographic information
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

    LOG.info("Found %s facets: %s", len(facets), facets)
    LOG.info("Processing %s time steps, saving every %s steps", len(size_df), n_steps)

    # Store results for analysis
    results_data = []

    # Process every n-th time step
    for i in range(0, len(size_df), n_steps):
        timestep = i
        row = size_df.iloc[i]

        # Extract energies (skip the 'time' column)
        energies = np.asarray(row.iloc[1:].dropna(), dtype=float)

        # Ensure we have matching facets and energies
        if len(energies) != len(facets):
            LOG.warning(
                "Mismatch at timestep %s: %s facets, %s energies",
                timestep,
                len(facets),
                len(energies),
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
                LOG.warning("Wulff shape file already exists, skipping: %s", output_filepath)
            else:
                wulff_result.export(str(output_filepath))
                LOG.debug("Saved Wulff shape: %s", output_filepath)

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

            LOG.debug("Processed timestep %s: %s facets", timestep, len(facets_subset))

        # RuntimeError covers scipy's QhullError, raised for degenerate hulls
        except (ValueError, IndexError, OSError, RuntimeError) as e:
            LOG.error("Error processing timestep %s: %s", timestep, e)
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
        with open(detailed_json_path, "w", encoding="utf-8") as f:
            json.dump(results_data, f, indent=2, cls=NumpyEncoder)

        LOG.info("Saved summary: %s", results_csv_path)
        LOG.info("Saved detailed data: %s", detailed_json_path)

    return results_data


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder for numpy arrays."""

    def default(self, o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        return super().default(o)


def process_multiple_size_files(
    size_files: list[Path],
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
    calculator = CrystalWulff.from_json(crystal_json_path)

    for size_file in size_files:
        if not isinstance(size_file, Path):
            size_file = Path(size_file)

        LOG.info("Processing size file: %s", size_file)

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
                    LOG.info("Saved combined Wulff processing results: %s", csv_output_path)
                else:
                    LOG.warning("No valid results to save to CSV")

            else:
                LOG.warning("No results generated to save")

        except (ValueError, IndexError, OSError, RuntimeError) as e:
            LOG.error("Failed to process %s: %s", size_file, e)
            continue

    return all_results


def build_parser() -> argparse.ArgumentParser:
    """Command line interface for the Wulff shape tools."""
    parser = argparse.ArgumentParser(
        prog="surfaces",
        description="Build Wulff shapes from a surface energies JSON file.",
    )
    parser.add_argument("json", type=Path, help="Surface energies JSON file")

    sub = parser.add_subparsers(dest="cmd", required=True)

    p_info = sub.add_parser("info", help="Print the crystal and surface cut summary")
    p_info.add_argument(
        "--context",
        default="solvated",
        choices=["vacuum", "solvated"],
        help="Energies used to pick the minimum energy cut (default: solvated)",
    )

    p_mesh = sub.add_parser("mesh", help="Export the Wulff shape as a mesh file")
    p_mesh.add_argument("-o", "--savedir", type=Path, default=Path("."), help="Output directory")
    p_mesh.add_argument("--name", default="wulff", help="Base name of the output file")
    p_mesh.add_argument("--suffix", default=".stl", help="Mesh format suffix (default: .stl)")

    p_xyz = sub.add_parser("xyz", help="Export the Wulff shape as an xyz point cloud")
    p_xyz.add_argument("-o", "--savedir", type=Path, default=Path("."), help="Output directory")
    p_xyz.add_argument("--name", default="wulff", help="Base name of the output file")
    p_xyz.add_argument("--suffix", default="txt", help="File extension (default: txt)")
    p_xyz.add_argument("--index", type=int, default=0, help="Index appended to the filename")
    p_xyz.add_argument(
        "--l-max",
        type=int,
        default=20,
        help="Maximum degree of the spherical harmonic expansion",
    )

    for p in (p_mesh, p_xyz):
        p.add_argument(
            "--solvation",
            default="solvated",
            choices=["vacuum", "solvated"],
            help="Which stored energies to use (default: solvated)",
        )
        p.add_argument(
            "--reduced",
            action="store_true",
            help="Keep only the lowest energy cut of each unique Miller index",
        )

    p_size = sub.add_parser("size", help="Generate Wulff shapes from size file(s)")
    p_size.add_argument("size_files", nargs="+", type=Path, help="Size CSV file(s)")
    p_size.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        type=Path,
        default=Path("."),
        help="Directory the combined results CSV is written to",
    )
    p_size.add_argument(
        "-n",
        "--n-steps",
        dest="n_steps",
        type=int,
        default=10,
        help="Generate a Wulff shape every n time steps (default: 10)",
    )
    p_size.add_argument(
        "--format",
        dest="file_format",
        default="stl",
        help="Mesh format of the generated shapes (default: stl)",
    )
    p_size.add_argument(
        "--no-expand-symmetry",
        dest="expand_symmetry",
        action="store_false",
        help="Do not expand facets using the symmetry operations",
    )
    p_size.add_argument(
        "--no-reduce-facets",
        dest="reduce_facets",
        action="store_false",
        help="Do not reduce facets to unique crystallographic directions",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "size":
        process_multiple_size_files(
            args.size_files,
            str(args.json),
            str(args.output_dir),
            n_steps=args.n_steps,
            file_format=args.file_format,
            expand_symmetry=args.expand_symmetry,
            reduce_facets=args.reduce_facets,
        )
        return 0

    calculator = CrystalWulff.from_json(args.json)

    if args.cmd == "info":
        info = calculator.get_crystal_info()
        print(f"Symmetry operations: {info['n_symmetry_operations']}")
        if calculator.hkl is None:
            print(f"No surface cut data in {args.json}")
            return 0

        print(f"Surface cuts: {calculator.n_surfaces}")
        calculator.log_surfaces(msg=str(args.json), show_offset=True, console=True)

        reduced = calculator.reduced_surface_energies(context=args.context)
        print(f"\nUnique Miller indices ({args.context}):")
        for hkl, vacuum, solvated in zip(
            reduced["hkl"], reduced["vacuum"], reduced["solvated"], strict=True
        ):
            print(f"{np.array2string(hkl):>12}  vacuum {vacuum:10.4f}  solvated {solvated:10.4f}")
        return 0

    if args.reduced:
        calculator = calculator.reduced_copy(context=args.solvation)

    if args.cmd == "mesh":
        mesh = calculator.to_mesh(
            name=args.name,
            savedir=args.savedir,
            suffix=args.suffix,
            solvation=args.solvation,
        )
        print(f"Wrote {args.savedir / (args.name + args.suffix)} (volume {mesh.volume:.4f})")
    elif args.cmd == "xyz":
        calculator.to_xyz(
            name=args.name,
            savedir=args.savedir,
            suffix=args.suffix,
            index=args.index,
            l_max=args.l_max,
            solvation=args.solvation,
        )
        print(f"Wrote {args.savedir / f'{args.name}_{args.index}.{args.suffix}'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
