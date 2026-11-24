"""Describe and draw the dimers behind each group of interaction energies.

A net file names the dimer behind every interaction line (e.g. ``dimer_6``),
and the energy calculation writes those dimers as XYZ files to a folder whose
name ends in ``dimers`` (e.g. ``ce-b3lyp_dimers``). :func:`export_key_dimers`
copies the dimer of each energy group to ``RESULTS/key_dimers`` named after its
group and dimer (``Group_4_dimer_0.xyz``), writes a lookup table of the contacts
holding each dimer together, and draws every group's dimer in one figure. Given
the crystal's CIF it also says which way each interaction runs in the crystal
and what it builds (see :mod:`cgpytools.analysis.motifs`).

Contacts are found from the geometry alone, with the usual crystallographic
cut-offs, so they are a quick description rather than an energy decomposition:

- hydrogen bonds: D-H...A with D, A in N, O, F; H...A <= 2.6 A, angle >= 120 deg
- weak hydrogen bonds: C-H...A with A in N, O; H...A <= 2.7 A, angle >= 130 deg
- halogen bonds: C-X...A with X in Cl, Br, I; X...A <= sum of vdW radii, angle >= 150 deg
- pi stacks: planar 5/6-rings, centroids <= 4.4 A apart, planes within 25 deg
  (offset if the centroids slip more than 1.6 A sideways)
- edge-to-face / C-H...pi: H <= 3.0 A from a ring centroid of the other molecule
"""

from __future__ import annotations

import logging
import re
import shutil
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cgpytools.analysis.motifs import CrystalMotifs, Pair

LOG = logging.getLogger("CG-ANALYSE")

COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Si": 1.11,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
}
VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80,
    "Cl": 1.75, "Br": 1.85, "I": 1.98,
}
ELEMENT_COLOURS = {
    "H": "#FFFFFF", "C": "#909090", "N": "#3050F8", "O": "#FF0D0D", "F": "#90E050",
    "P": "#FF8000", "S": "#E6C200", "Cl": "#1FF01F", "Br": "#A62929", "I": "#940094",
}
CONTACT_STYLES = {
    "hbond": {"color": "#E07B00", "linestyle": (0, (4, 2)), "label": "H-bond"},
    "weak": {"color": "#C9A227", "linestyle": (0, (1, 2)), "label": "weak C-H...O/N"},
    "halogen": {"color": "#1E9E50", "linestyle": (0, (4, 2)), "label": "halogen bond"},
    "pi": {"color": "#7B3FBF", "linestyle": (0, (2, 2)), "label": "pi stack"},
    "chpi": {"color": "#3F8FBF", "linestyle": (0, (1, 2)), "label": "C-H...pi"},
}
HBOND_DONORS = {"N", "O", "F"}
HBOND_ACCEPTORS = {"N", "O", "F"}
HALOGENS = {"Cl", "Br", "I"}
RING_ELEMENTS = {"C", "N", "O", "S"}
MULTIPLICITY = {1: "Single", 2: "Double", 3: "Triple", 4: "Quadruple"}


@dataclass
class Contact:
    kind: str  # a key of CONTACT_STYLES
    label: str  # e.g. "N-H...O"
    distance: float
    points: tuple[np.ndarray, np.ndarray]  # ends of the line to draw


@dataclass
class Dimer:
    elements: list[str]
    coords: np.ndarray
    molecule: np.ndarray  # 0 or 1 for each atom
    bonds: list[tuple[int, int]]
    contacts: list[Contact] = field(default_factory=list)
    closest: float = np.nan
    closest_pair: str = ""

    @property
    def description(self) -> str:
        """A short description of the contacts, e.g. ``Double H-bond (N-H...O x2)``."""
        counts = {kind: Counter() for kind in CONTACT_STYLES}
        for contact in self.contacts:
            counts[contact.kind][contact.label] += 1

        def listed(counter: Counter) -> str:
            return ", ".join(f"{label} x{n}" if n > 1 else label for label, n in counter.items())

        parts = []
        n_hbonds = sum(counts["hbond"].values())
        if n_hbonds:
            word = MULTIPLICITY.get(n_hbonds, f"{n_hbonds}x")
            parts.append(f"{word} H-bond ({listed(counts['hbond'])})")
        if counts["halogen"]:
            parts.append(f"halogen bond ({listed(counts['halogen'])})")
        stacks = [c for c in self.contacts if c.kind == "pi"]
        if stacks:
            kinds = Counter(c.label for c in stacks)
            nearest = min(c.distance for c in stacks)
            parts.append(f"{listed(kinds)} ({nearest:.2f} A)")
        if counts["chpi"]:
            parts.append(f"edge-to-face C-H...pi x{sum(counts['chpi'].values())}")
        if counts["weak"]:
            parts.append(f"weak {listed(counts['weak'])}")
        if not parts:
            return f"Dispersion only (closest {self.closest_pair} {self.closest:.2f} A)"
        text = "; ".join(parts)
        return text[0].upper() + text[1:]


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    """Elements and (n, 3) coordinates of the first frame of an XYZ file."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    n_atoms = int(lines[0].split()[0])
    elements, coords = [], []
    for line in lines[2 : 2 + n_atoms]:
        parts = line.split()
        symbol = "".join(ch for ch in parts[0] if ch.isalpha())
        elements.append(symbol[:1].upper() + symbol[1:].lower())
        coords.append([float(value) for value in parts[1:4]])
    return elements, np.array(coords)


def _find_bonds(elements: list[str], coords: np.ndarray) -> list[tuple[int, int]]:
    radii = np.array([COVALENT_RADII.get(el, 0.8) for el in elements])
    dist = np.linalg.norm(coords[:, None] - coords[None], axis=-1)
    limit = (radii[:, None] + radii[None]) * 1.15
    i, j = np.nonzero(np.triu((dist < limit) & (dist > 0.1), 1))
    return list(zip(i.tolist(), j.tolist()))


def _split_molecules(n_atoms: int, bonds: list[tuple[int, int]]) -> np.ndarray:
    """Label each atom 0 or 1 for the two molecules of the dimer.

    Bonded fragments stay together; fragments are assigned to the molecule
    whose half of the file they start in (dimer files list one molecule, then
    the other), which also handles salts and solvates.
    """
    parent = list(range(n_atoms))

    def root(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in bonds:
        parent[root(a)] = root(b)
    fragments: dict[int, list[int]] = {}
    for atom in range(n_atoms):
        fragments.setdefault(root(atom), []).append(atom)

    molecule = np.zeros(n_atoms, dtype=int)
    if len(fragments) == 2:
        second = max(fragments.values(), key=min)
        molecule[second] = 1
    else:
        for atoms in fragments.values():
            molecule[atoms] = int(min(atoms) >= n_atoms / 2)
    return molecule


def _rings(elements, coords, bonds, atoms: set[int]) -> list[tuple[np.ndarray, np.ndarray]]:
    """Centroid and unit normal of each planar 5- or 6-membered ring among ``atoms``."""
    neighbours = {a: set() for a in atoms}
    for a, b in bonds:
        if a in atoms and b in atoms and elements[a] != "H" and elements[b] != "H":
            neighbours[a].add(b)
            neighbours[b].add(a)

    found = set()
    def walk(path):
        if len(path) > 6:
            return
        for nxt in neighbours[path[-1]]:
            if nxt == path[0] and len(path) >= 5:
                found.add(frozenset(path))
            elif nxt not in path and nxt > path[0]:
                walk(path + [nxt])

    for start in atoms:
        if elements[start] in RING_ELEMENTS:
            walk([start])

    rings = []
    for ring in found:
        if len(ring) not in (5, 6) or any(elements[a] not in RING_ELEMENTS for a in ring):
            continue
        points = coords[list(ring)]
        centroid = points.mean(axis=0)
        _, singular, vt = np.linalg.svd(points - centroid)
        normal = vt[2]
        if np.abs((points - centroid) @ normal).max() < 0.15:  # planar, so likely aromatic
            rings.append((centroid, normal))
    return rings


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Angle a-b-c in degrees."""
    u, v = a - b, c - b
    cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(cos, -1, 1))))


def analyse_dimer(path: Path) -> Dimer:
    """Read a dimer XYZ file and find the contacts between its two molecules."""
    elements, coords = read_xyz(path)
    bonds = _find_bonds(elements, coords)
    molecule = _split_molecules(len(elements), bonds)
    dimer = Dimer(elements, coords, molecule, bonds)

    bonded_to = {a: [] for a in range(len(elements))}
    for a, b in bonds:
        bonded_to[a].append(b)
        bonded_to[b].append(a)

    halves = [set(np.flatnonzero(molecule == m).tolist()) for m in (0, 1)]
    dist = np.linalg.norm(coords[:, None] - coords[None], axis=-1)
    cross = dist[np.ix_(sorted(halves[0]), sorted(halves[1]))]
    if cross.size:
        i, j = np.unravel_index(np.argmin(cross), cross.shape)
        a, b = sorted(halves[0])[i], sorted(halves[1])[j]
        dimer.closest = float(cross[i, j])
        dimer.closest_pair = f"{elements[a]}...{elements[b]}"

    for here, there in ((halves[0], halves[1]), (halves[1], halves[0])):
        for h in here:
            if elements[h] != "H" or not bonded_to[h]:
                continue
            donor = bonded_to[h][0]
            for acceptor in there:
                el = elements[acceptor]
                d = dist[h, acceptor]
                if el not in HBOND_ACCEPTORS:
                    continue
                angle = _angle(coords[donor], coords[h], coords[acceptor])
                label = f"{elements[donor]}-H...{el}"
                if elements[donor] in HBOND_DONORS and d <= 2.6 and angle >= 120:
                    dimer.contacts.append(Contact("hbond", label, d, (coords[h], coords[acceptor])))
                elif elements[donor] == "C" and d <= 2.7 and angle >= 130:
                    dimer.contacts.append(Contact("weak", label, d, (coords[h], coords[acceptor])))

        for x in here:
            if elements[x] not in HALOGENS or not bonded_to[x]:
                continue
            carbon = bonded_to[x][0]
            for acceptor in there:
                el = elements[acceptor]
                if el not in {"N", "O", "S", *HALOGENS}:
                    continue
                d = dist[x, acceptor]
                limit = VDW_RADII.get(elements[x], 1.8) + VDW_RADII.get(el, 1.7)
                if d <= limit and _angle(coords[carbon], coords[x], coords[acceptor]) >= 150:
                    label = f"C-{elements[x]}...{el}"
                    ends = (coords[x], coords[acceptor])
                    dimer.contacts.append(Contact("halogen", label, d, ends))

    rings = [_rings(elements, coords, bonds, half) for half in halves]
    for (c1, n1), (c2, n2) in ((r1, r2) for r1 in rings[0] for r2 in rings[1]):
        offset = c2 - c1
        d = float(np.linalg.norm(offset))
        tilt = np.degrees(np.arccos(min(1.0, abs(float(np.dot(n1, n2))))))
        if d <= 4.4 and tilt <= 25:
            slip = float(np.sqrt(max(d**2 - np.dot(offset, n1) ** 2, 0)))
            label = "offset pi stack" if slip > 1.6 else "face-to-face pi stack"
            dimer.contacts.append(Contact("pi", label, d, (c1, c2)))

    for m, other in ((0, 1), (1, 0)):
        for centroid, _ in rings[other]:
            for h in halves[m]:
                if elements[h] == "H" and np.linalg.norm(coords[h] - centroid) <= 3.0:
                    d = float(np.linalg.norm(coords[h] - centroid))
                    dimer.contacts.append(Contact("chpi", "C-H...pi", d, (coords[h], centroid)))
    return dimer


def _view_basis(dimer: Dimer) -> np.ndarray:
    """Screen axes: x joins the two molecules, y spans most of the rest, z is depth."""
    centroids = [dimer.coords[dimer.molecule == m].mean(axis=0) for m in (0, 1)]
    x = centroids[1] - centroids[0]
    if np.linalg.norm(x) < 1e-6:
        x = np.array([1.0, 0.0, 0.0])
    x /= np.linalg.norm(x)
    centred = dimer.coords - dimer.coords.mean(axis=0)
    rest = centred - np.outer(centred @ x, x)
    _, _, vt = np.linalg.svd(rest)
    y = vt[0] - np.dot(vt[0], x) * x
    y /= np.linalg.norm(y)
    return np.vstack([x, y, np.cross(x, y)])


AXIS_COLOURS = {"a": "#D62728", "b": "#2CA02C", "c": "#1F77B4"}


def _draw_lattice_axes(ax, basis: np.ndarray, lattice: np.ndarray):
    """Arrows for the a, b and c axes as seen in this view, in the panel's corner.

    Arrow length shows how much of each axis lies in the page. An axis
    (nearly) along the line of sight is listed beside the origin instead, as
    ⊙ when it points out of the page and ⊗ when it points in.
    """
    origin = np.array([0.1, 0.12])
    length = 0.1
    end_on, labels = [], []
    for name, axis in zip("abc", lattice):
        unit = axis / np.linalg.norm(axis)
        dx, dy, dz = basis @ unit
        colour = AXIS_COLOURS[name]
        in_page = np.hypot(dx, dy)
        if in_page < 0.25:
            end_on.append((name, colour, "⊙" if dz > 0 else "⊗"))
            continue
        tip = origin + length * np.array([dx, dy])
        ax.annotate("", xy=tip, xytext=origin, xycoords="axes fraction",
                    arrowprops={"arrowstyle": "-|>", "color": colour, "lw": 1.3,
                                "shrinkA": 0, "shrinkB": 0, "mutation_scale": 9})
        label = tip + 0.045 * np.array([dx, dy]) / in_page
        # Two axes can share a projection (looking down [1-10], a and b do):
        # set their labels side by side
        while any(np.linalg.norm(label - other) < 0.04 for other in labels):
            label = label + 0.04 * np.array([dy, -dx]) / in_page
        labels.append(label)
        ax.text(*label, name, transform=ax.transAxes, color=colour, fontsize=9,
                fontstyle="italic", fontweight="bold", ha="center", va="center", zorder=6)
    ax.plot(*origin, "o", color="#404040", markersize=2.5, transform=ax.transAxes, zorder=6)
    for i, (name, colour, symbol) in enumerate(end_on):
        ax.text(origin[0] + 0.035, origin[1] - 0.045 - 0.055 * i, f"{symbol} {name}",
                transform=ax.transAxes, color=colour, fontsize=9, fontstyle="italic",
                fontweight="bold", ha="left", va="center", zorder=6)


def draw_dimer(ax, dimer: Dimer, title: str = "", subtitle: str = "",
               lattice: np.ndarray | None = None):
    """Ball-and-stick picture of a dimer with its intermolecular contacts marked.

    ``lattice`` (rows a, b, c in the dimer file's Cartesian frame) adds
    arrows showing the crystal axes in this view.
    """
    basis = _view_basis(dimer)
    centre = dimer.coords.mean(axis=0)
    xyz = (dimer.coords - centre) @ basis.T
    project = lambda p: (np.asarray(p) - centre) @ basis.T  # noqa: E731

    bond_colours = ("#404040", "#2F5F8F")
    for a, b in dimer.bonds:
        colour = bond_colours[dimer.molecule[a]]
        ax.plot(*xyz[[a, b], :2].T, color=colour, linewidth=2.2, zorder=2, solid_capstyle="round")

    for contact in dimer.contacts:
        style = CONTACT_STYLES[contact.kind]
        ends = np.array([project(p) for p in contact.points])
        ax.plot(*ends[:, :2].T, color=style["color"], linestyle=style["linestyle"],
                linewidth=1.6, zorder=3)

    order = np.argsort(xyz[:, 2])  # far atoms first
    for atom in order:
        el = dimer.elements[atom]
        size = 22 if el == "H" else 70
        ax.scatter(*xyz[atom, :2], s=size, color=ELEMENT_COLOURS.get(el, "#FF69B4"),
                   edgecolor=bond_colours[dimer.molecule[atom]], linewidth=0.8, zorder=4)

    ax.set_aspect("equal")
    ax.margins(0.08)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if lattice is not None:
        _draw_lattice_axes(ax, basis, lattice)
    if title:
        ax.set_title(title, fontsize=11, fontweight="bold", pad=4)
    if subtitle:
        ax.text(0.5, -0.02, subtitle, transform=ax.transAxes, ha="center", va="top",
                fontsize=8.5, wrap=True)


def find_dimers_folder(input_dir: Path | None) -> Path | None:
    """The folder of dimer XYZ files in ``input_dir``: one whose name ends in 'dimers'."""
    if input_dir is None or not Path(input_dir).is_dir():
        return None

    def is_source(p: Path) -> bool:
        # key_dimers is this module's own output, so never a source
        name = p.name.lower()
        return p.is_dir() and name.endswith("dimers") and name != "key_dimers"

    candidates = [p for p in Path(input_dir).iterdir() if is_source(p)]
    if not candidates:
        candidates = [p for p in Path(input_dir).rglob("*") if is_source(p)]
    return min(candidates, key=lambda p: len(p.parts), default=None)


def export_key_dimers(
    groups: dict[str, list[str]],
    dimers_folder: Path,
    out_folder: Path,
    energies: pd.DataFrame | None = None,
    multiplicity: dict[str, int] | None = None,
    distances: dict[str, float] | None = None,
    pairs: dict[str, list[Pair | None]] | None = None,
    crystal: CrystalMotifs | None = None,
) -> pd.DataFrame:
    """Copy, describe and draw the dimer behind each energy group.

    ``groups`` maps a group name (``Group_2``) to the dimer names of its net
    file lines (``["dimer_4", "dimer_4"]``). ``energies`` holds the group
    energy columns, summarised across the batch in the lookup table.
    With ``crystal`` and ``pairs`` (the parsed molecule pair of each of those
    lines) each dimer is also placed in the crystal: the lattice direction of
    its interaction vector and the motif the interaction builds (chain, layer...).
    Returns the lookup table, also saved as ``key_dimers.csv``.
    """
    out_folder.mkdir(parents=True, exist_ok=True)
    files = {p.stem.lower(): p for p in dimers_folder.glob("*.xyz")}

    rows, drawn = [], []
    for group, dimer_names in groups.items():
        names = list(dict.fromkeys(n for n in dimer_names if n))
        for name in names:
            source = files.get(name.lower())
            if source is None:
                LOG.warning("No %s.xyz in %s for %s", name, dimers_folder, group)
                continue
            target = out_folder / f"{group}_{name}.xyz"
            shutil.copy2(source, target)

            row = {"group": group, "dimer": name, "file": target.name}
            if multiplicity:
                row["lines_in_net"] = multiplicity.get(group)
            if distances:
                row["R (A)"] = distances.get(group)
            try:
                dimer = analyse_dimer(source)
                row["contacts"] = dimer.description
                row["closest contact (A)"] = round(dimer.closest, 2)
                drawn.append((group, name, dimer))
            except Exception:
                LOG.exception("Could not read dimer %s", source)
                row["contacts"] = "unreadable"
                dimer = None
            if crystal is not None and pairs:
                lines = [
                    p for p, n in zip(pairs.get(group, []), dimer_names) if p is not None
                    and n == name
                ]
                if lines:
                    _place_in_crystal(row, crystal, lines, dimer)
            if energies is not None and group in energies.columns:
                values = pd.to_numeric(energies[group], errors="coerce")
                row["energy mean"] = round(values.mean(), 3)
                row["energy min"] = round(values.min(), 3)
                row["energy max"] = round(values.max(), 3)
            rows.append(row)

    table = pd.DataFrame(rows)
    table = table.drop(columns=[c for c in table.columns if c.startswith("_")])
    if table.empty:
        LOG.warning("No key dimers found in %s", dimers_folder)
        return table
    table.to_csv(out_folder / "key_dimers.csv", index=False)

    LOG.info("Key dimers (%s):", out_folder)
    for row in rows:
        LOG.info("  %-10s %-10s %s", row["group"], row["dimer"], row.get("contacts", ""))
        if "motif" in row:
            LOG.info("  %21s vector %s; %s", "", row["direction"], row["motif"])

    if drawn:
        _draw_key_dimers(drawn, rows, out_folder / "key_dimers.png", crystal)
    return table


def _place_in_crystal(row: dict, crystal: CrystalMotifs, lines: list[Pair], dimer: Dimer | None):
    """Add the interaction direction and motif of a dimer to its table row.

    The direction is that of the drawn dimer when its file is in the CIF's
    frame (so it matches the picture's a/b/c arrows), else that of the
    first net file line naming it; both are copies of the same interaction.
    """
    try:
        vector = None
        if dimer is not None:
            vector = crystal.dimer_vector(dimer.coords, dimer.molecule)
        row["_in_frame"] = vector is not None
        if vector is None:
            vector = crystal.vector(lines[0])
        direction = crystal.direction(vector)
        motif = crystal.motif(lines)
    except Exception:
        LOG.exception("Could not place %s in the crystal", row["dimer"])
        return
    row["vector (fractional)"] = "(" + ", ".join(f"{v:.3f}" for v in direction.vector) + ")"
    row["direction"] = direction.text
    row["angle to a/b/c (deg)"] = "/".join(f"{a:.0f}" for a in direction.axis_angles)
    row["motif"] = motif.text
    row["_direction_math"] = direction.math_text
    row["_motif_math"] = motif.math_text


def _pretty(text: str) -> str:
    """Typeset a plain-ASCII contact description for the figure."""
    text = text.replace("...", "···").replace(" A)", " Å)")
    text = re.sub(r"\bpi\b", "π", text)
    return re.sub(r" x(\d)", r" ×\1", text)


def _draw_key_dimers(drawn, rows, path: Path, crystal: CrystalMotifs | None = None):
    n = len(drawn)
    n_cols = min(4, n)
    n_rows = int(np.ceil(n / n_cols))
    placed = any("_motif_math" in row for row in rows)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(4.2 * n_cols, (4.8 if placed else 4.4) * n_rows + 0.6),
        facecolor="white", squeeze=False,
    )
    by_key = {(row["group"], row["dimer"]): row for row in rows}
    for ax, (group, name, dimer) in zip(axes.flat, drawn):
        row = by_key[(group, name)]
        title = f"{group.replace('_', ' ')}  ({name})"
        lines = [_pretty(dimer.description)]
        if "energy mean" in row:
            lines.append(f"E = {row['energy mean']:.2f} (range {row['energy min']:.2f} to "
                         f"{row['energy max']:.2f}) kJ/mol")
        if row.get("R (A)") is not None:
            lines.append(f"R = {row['R (A)']:.2f} Å")
        if "_motif_math" in row:
            lines.append(f"Vector {row['_direction_math']}; {row['_motif_math']}")
        lattice = crystal.direct if crystal is not None and row.get("_in_frame") else None
        draw_dimer(ax, dimer, title=title, subtitle="\n".join(lines), lattice=lattice)
    for ax in axes.flat[n:]:
        ax.set_visible(False)

    # One scale for every panel, so the dimers' sizes and separations compare directly
    half = max(
        max(np.ptp(ax.get_xlim()), np.ptp(ax.get_ylim())) for ax in axes.flat[:n]
    ) / 2
    for ax in axes.flat[:n]:
        cx, cy = np.mean(ax.get_xlim()), np.mean(ax.get_ylim())
        ax.set_xlim(cx - half, cx + half)
        ax.set_ylim(cy - half, cy + half)

    used = {c.kind for _, _, d in drawn for c in d.contacts}
    handles = [
        plt.Line2D(
            [], [], color=s["color"], linestyle=s["linestyle"], linewidth=1.6,
            label=_pretty(s["label"]),
        )
        for kind, s in CONTACT_STYLES.items()
        if kind in used
    ]
    if handles:
        fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False)
    fig.suptitle("Key dimers by energy group", fontsize=14, fontweight="bold")
    # Room below the last row for its subtitle lines and the legend
    bottom = (0.4 + (0.25 if placed else 0)) / fig.get_figheight() if handles else 0
    fig.tight_layout(rect=(0, bottom, 1, 0.97), h_pad=4.5)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    LOG.info("Saved %s", path)


__all__ = ["analyse_dimer", "draw_dimer", "export_key_dimers", "find_dimers_folder"]
