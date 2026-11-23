"""Draw crystal shapes as small shaded 2D pictures for annotating plots.

A shape (a CrystoGen particle cloud or a Wulff mesh) is drawn as the convex hull
of its coordinates, seen from a fixed oblique angle relative to the crystal's
own principal axes, so pictures of different crystals are directly comparable:
the longest axis runs left-right, the medium axis up-down, and the shortest
axis points roughly towards the viewer, tilted so the thickness is visible.
"""

from pathlib import Path

import numpy as np
import trimesh
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import to_rgb
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from cgpytools.io.crystal import CrystalShape

# Viewing angle applied after principal-axis alignment, in degrees
VIEW_TILT = 30.0
VIEW_TURN = -35.0
# Light direction in view space (towards the viewer, from the upper left)
LIGHT = np.array([-0.45, 0.6, 0.65]) / np.linalg.norm([-0.45, 0.6, 0.65])
# Neighbouring hull faces closer than this angle are drawn as one facet
CREASE_ANGLE = 10.0


def load_shape_points(path: Path, frame: int = -1) -> np.ndarray:
    """Return the (n, 3) coordinates of one frame of a shape file.

    Accepts every format :meth:`CrystalShape.from_file` reads (XYZ, TXT, STL,
    PLY, GLB). Mesh formats hold a single frame.
    """
    crystal = CrystalShape.from_file(Path(path), normalise=False)
    return np.asarray(crystal.frames[frame].coords, dtype=float)


def _rotation(axis: str, degrees: float) -> np.ndarray:
    c, s = np.cos(np.radians(degrees)), np.sin(np.radians(degrees))
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def view_coordinates(points: np.ndarray) -> np.ndarray:
    """Centre, align to principal axes, scale to unit size and apply the view angle."""
    centred = points - points.mean(axis=0)
    _, _, vh = np.linalg.svd(centred, full_matrices=False)
    if np.linalg.det(vh) < 0:
        # Keep the basis right-handed so the crystal is not drawn mirrored
        vh[2] *= -1
    aligned = centred @ vh.T
    aligned /= np.linalg.norm(aligned, axis=1).max() or 1.0
    view = _rotation("y", VIEW_TURN) @ _rotation("x", VIEW_TILT)
    return aligned @ view.T


def _shaded_colors(base: str, intensity: np.ndarray) -> np.ndarray:
    """Blend ``base`` from a dark to a light tint by ``intensity`` in [0, 1]."""
    rgb = np.array(to_rgb(base))
    dark = rgb * 0.55
    light = rgb + (1 - rgb) * 0.6
    return dark + (light - dark) * intensity[:, None]


def draw_crystal(
    ax,
    points: np.ndarray,
    color: str = "#2E86AB",
    edge_color: str = "#263238",
    crease_angle: float = CREASE_ANGLE,
) -> bool:
    """Draw the convex hull of ``points`` into ``ax`` as a shaded solid.

    Only the faces turned towards the viewer are drawn (the hull is convex, so
    they never overlap). Facet edges and the outline are drawn as lines;
    near-coplanar hull triangles are merged into one flat-shaded facet, so a
    rough particle-cloud surface still reads as clean crystal faces.

    Returns False, leaving a short note in ``ax``, if no 3D hull can be built.
    """
    ax.set_axis_off()
    ax.set_aspect("equal")

    try:
        mesh = trimesh.convex.convex_hull(view_coordinates(points))
    except Exception:  # noqa: BLE001 - qhull raises several types for degenerate input
        mesh = None
    if mesh is None or len(mesh.faces) < 4:
        ax.text(0.5, 0.5, "no 3D shape", ha="center", va="center", transform=ax.transAxes)
        return False

    normals = mesh.face_normals
    pairs = mesh.face_adjacency
    smooth = mesh.face_adjacency_angles < np.radians(crease_angle)

    # Group faces joined by smooth edges into facets with one area-weighted normal
    n_faces = len(mesh.faces)
    graph = coo_matrix(
        (np.ones(smooth.sum()), (pairs[smooth, 0], pairs[smooth, 1])), shape=(n_faces, n_faces)
    )
    _, facet = connected_components(graph, directed=False)
    facet_normals = np.zeros((facet.max() + 1, 3))
    np.add.at(facet_normals, facet, normals * mesh.area_faces[:, None])
    facet_normals /= np.linalg.norm(facet_normals, axis=1, keepdims=True)
    shading = np.clip(facet_normals[facet] @ LIGHT, 0, 1)

    front = normals[:, 2] > 0
    colors = _shaded_colors(color, shading[front])
    polygons = mesh.vertices[mesh.faces[front]][:, :, :2]
    # Edges in the face colour hide the antialiasing seams between triangles of one facet
    ax.add_collection(
        PolyCollection(
            polygons, facecolors=colors, edgecolors=colors, linewidths=1.0, joinstyle="round"
        )
    )

    front_a, front_b = front[pairs[:, 0]], front[pairs[:, 1]]
    visible_crease = ~smooth & (front_a | front_b)
    silhouette = front_a != front_b
    edges = mesh.vertices[mesh.face_adjacency_edges[visible_crease | silhouette]][:, :, :2]
    ax.add_collection(
        LineCollection(edges, colors=edge_color, linewidths=0.9, capstyle="round")
    )

    lo, hi = mesh.vertices[:, :2].min(axis=0), mesh.vertices[:, :2].max(axis=0)
    pad = 0.05 * (hi - lo).max()
    ax.set_xlim(lo[0] - pad, hi[0] + pad)
    ax.set_ylim(lo[1] - pad, hi[1] + pad)
    return True
