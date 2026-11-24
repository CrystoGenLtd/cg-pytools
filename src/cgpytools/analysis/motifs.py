"""Place the dimers of a net file in the crystal: their direction and the motif they build.

Every interaction line of a net file names its two molecules by symmetry
operation, e.g. ``[1A(y,-x,-z+1)-1A(x,y,z+1)]``. With the crystal's CIF this
gives, for each energy group:

- the **interaction vector** between the molecule centroids, as a fractional
  vector and the nearest low-index lattice direction ``[uvw]``;
- the **motif** that interaction builds on its own when every symmetry copy of
  it is applied: a finite dimer, ring or cluster, a chain along ``[uvw]``, a
  layer parallel to ``(hkl)``, or a 3D network.

The motif is found from the periodic graph of the molecules joined by that one
interaction: a cycle in the unit cell graph that closes through a lattice
translation is a chain along that translation, and the rank of all such
translations gives the dimensionality (0 finite, 1 chain, 2 layer, 3 network).
"""

from __future__ import annotations

import itertools
import logging
import re
from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from math import gcd
from pathlib import Path

import numpy as np

LOG = logging.getLogger("CG-ANALYSE")

PAIR_RE = re.compile(r"\[(\d[A-Za-z])(?:\(([^)]*)\))?-(\d[A-Za-z])(?:\(([^)]*)\))?\]")
R_TOLERANCE = 0.05  # Å between a net file's R and the centroid distance it implies
FRAME_TOLERANCE = 0.15  # Å between a dimer file's centroid and a crystal site
MAX_INDEX = 3  # largest index tried when naming a direction [uvw]
ALONG = 1.0  # degrees within which a vector is called "along" a direction
NEAR = 5.0  # degrees within which it is called "~" one, else given by its angles to a, b, c


@dataclass(frozen=True)
class Op:
    """A symmetry operation x' = R x + t that keeps its whole translation (z+1 stays +1)."""

    rotation: np.ndarray
    translation: np.ndarray

    @classmethod
    def parse(cls, text: str | None) -> Op:
        rotation, translation = np.zeros((3, 3)), np.zeros(3)
        if not text:
            return cls(np.eye(3), translation)
        rows = text.lower().replace(" ", "").split(",")
        if len(rows) != 3:
            raise ValueError(f"Not a symmetry operation: {text!r}")
        for i, row in enumerate(rows):
            for term in re.findall(r"[+-]?[^+-]+", row):
                sign = -1.0 if term.startswith("-") else 1.0
                term = term.lstrip("+-")
                if term[-1] in "xyz":
                    coefficient = term[:-1].rstrip("*")
                    rotation[i, "xyz".index(term[-1])] += sign * float(Fraction(coefficient or 1))
                else:
                    translation[i] += sign * float(Fraction(term))
        return cls(rotation, translation)

    def __call__(self, point: np.ndarray) -> np.ndarray:
        return self.rotation @ point + self.translation


@dataclass(frozen=True)
class Pair:
    """The two molecules of a net file line, each a (label, operation)."""

    first: str
    first_op: Op
    second: str
    second_op: Op
    r: float

    @classmethod
    def parse(cls, info: str, r: float) -> Pair | None:
        match = PAIR_RE.search(info or "")
        if match is None:
            return None
        first, op1, second, op2 = match.groups()
        return cls(first, Op.parse(op1), second, Op.parse(op2), float(r))


@dataclass
class Motif:
    rank: int  # 0 finite, 1 chain, 2 layer, 3 network
    n_molecules: int = 0  # molecules in a finite motif
    ring: bool = False  # finite motif in which every molecule has two partners
    index: tuple[int, int, int] | None = None  # chain [uvw] or layer (hkl)
    repeat: float | None = None  # chain repeat in Å

    @property
    def text(self) -> str:
        if self.rank == 0:
            if self.n_molecules <= 2:
                return "discrete dimer"
            if self.ring:
                return f"closed {self.n_molecules}-membered ring"
            return f"finite {self.n_molecules}-molecule cluster"
        if self.rank == 1:
            return f"chain along {bracket(self.index)} ({self.repeat:.2f} A repeat)"
        if self.rank == 2:
            return f"layer parallel to {bracket(self.index, '()')}"
        return "3D network"

    @property
    def math_text(self) -> str:
        """``text`` with indices typeset for matplotlib (bars over negatives)."""
        if self.rank == 1:
            return f"chain along {bracket_math(self.index)} ({self.repeat:.2f} Å repeat)"
        if self.rank == 2:
            return f"layer ∥ {bracket_math(self.index, '()')}"
        return self.text


def bracket(index, brackets: str = "[]") -> str:
    """``[1 -1 0]``: an index as plain text."""
    return brackets[0] + " ".join(str(int(i)) for i in index) + brackets[1]


def bracket_math(index, brackets: str = "[]") -> str:
    """``$[1\\bar{1}0]$``: an index as matplotlib mathtext, negatives overbarred."""
    open_, close = (r"\left[", r"\right]") if brackets == "[]" else (r"\left(", r"\right)")
    digits = "".join(rf"\bar{{{-int(i)}}}" if i < 0 else str(int(i)) for i in index)
    return f"${open_}{digits}{close}$"


def _primitive(vector) -> tuple[int, int, int]:
    """An integer vector divided by its gcd, first non-zero entry positive."""
    vector = [round(float(v)) for v in vector]
    divisor = reduce(gcd, (abs(v) for v in vector)) or 1
    vector = [v // divisor for v in vector]
    first = next((v for v in vector if v), 1)
    return tuple(v if first > 0 else -v for v in vector)


@dataclass
class Direction:
    vector: np.ndarray  # fractional
    index: tuple[int, int, int]  # nearest low-index direction
    deviation: float  # degrees between vector and index
    axis_angles: tuple[float, float, float]  # degrees to a, b, c (0-90)

    @property
    def _angles(self) -> str:
        return "/".join(f"{a:.0f}" for a in self.axis_angles)

    @property
    def text(self) -> str:
        if self.deviation < ALONG:
            return f"along {bracket(self.index)}"
        if self.deviation < NEAR:
            return f"~{bracket(self.index)} ({self.deviation:.1f} deg off)"
        return f"at {self._angles} deg to a/b/c"

    @property
    def math_text(self) -> str:
        if self.deviation < ALONG:
            return f"along {bracket_math(self.index)}"
        if self.deviation < NEAR:
            return f"≈ {bracket_math(self.index)} ({self.deviation:.1f}° off)"
        return f"at {self._angles}° to $a/b/c$"


class CrystalMotifs:
    """The molecules of a net file placed in the crystal from its CIF."""

    def __init__(self, cif: Path, pairs: list[Pair]):
        from chmpy import Crystal

        crystal = Crystal.load(str(cif))
        self.cif = Path(cif)
        self.direct = np.asarray(crystal.unit_cell.direct, dtype=float)  # rows a, b, c
        self.inverse = np.linalg.inv(self.direct)
        ops = list(crystal.space_group.symmetry_operations)
        ops.sort(key=lambda s: not s.is_identity())
        self.symops = [Op(np.asarray(s.rotation, float), np.asarray(s.translation, float))
                       for s in ops]
        molecules = crystal.symmetry_unique_molecules()
        self._unique = [crystal.to_fractional(m.centroid.reshape(1, 3))[0] for m in molecules]
        self.centroids = self._fit_centroids([p for p in pairs if p is not None])
        self._sites = [
            (label, s(c)) for label, c in self.centroids.items() for s in self.symops
        ]

    # -- geometry -----------------------------------------------------------

    def cartesian(self, fractional: np.ndarray) -> np.ndarray:
        return np.asarray(fractional) @ self.direct

    def fractional(self, cartesian: np.ndarray) -> np.ndarray:
        return np.asarray(cartesian) @ self.inverse

    def _length(self, fractional: np.ndarray) -> float:
        return float(np.linalg.norm(self.cartesian(fractional)))

    def _mismatch(self, pair: Pair, centroids: dict[str, np.ndarray]) -> float:
        vector = pair.second_op(centroids[pair.second]) - pair.first_op(centroids[pair.first])
        return abs(self._length(vector) - pair.r)

    def _fit_centroids(self, pairs: list[Pair]) -> dict[str, np.ndarray]:
        """The centroid of each molecule label the net file's operations act on.

        The net file's ``1A`` is one copy of a symmetry-unique molecule; the
        copy is the first (identity, no shift preferred) that reproduces every
        R of the net file it is in.
        """
        if not pairs:
            raise ValueError("no interaction lines with symmetry operations")
        labels = sorted({p.first for p in pairs} | {p.second for p in pairs})
        candidates = [
            s(c) + np.array(shift)
            for c in self._unique
            for s in self.symops
            for shift in sorted(itertools.product((0, -1, 1), repeat=3), key=np.count_nonzero)
        ]
        fitted: dict[str, np.ndarray] = {}
        remaining = list(labels)
        while remaining:
            progress = False
            for label in list(remaining):
                lines = [
                    p for p in pairs
                    if label in (p.first, p.second)
                    and {p.first, p.second} <= set(fitted) | {label}
                ]
                if not lines:
                    continue
                for candidate in candidates:
                    trial = {**fitted, label: candidate}
                    if all(self._mismatch(p, trial) < R_TOLERANCE for p in lines):
                        fitted[label] = candidate
                        remaining.remove(label)
                        progress = True
                        break
                else:
                    raise ValueError(
                        f"no molecule in the CIF reproduces the net file's R values for {label}"
                    )
            if not progress:
                raise ValueError(f"could not place molecules {remaining}")
        LOG.debug("Net file molecules placed at %s", fitted)
        return fitted

    # -- directions -----------------------------------------------------------

    def vector(self, pair: Pair) -> np.ndarray:
        """Fractional vector from the first molecule's centroid to the second's."""
        return pair.second_op(self.centroids[pair.second]) - pair.first_op(
            self.centroids[pair.first]
        )

    def direction(self, vector: np.ndarray) -> Direction:
        """Name a fractional vector by its nearest low-index lattice direction."""
        cart = self.cartesian(vector)
        unit = cart / np.linalg.norm(cart)
        best, best_angle = (0, 0, 1), 180.0
        for index in itertools.product(range(-MAX_INDEX, MAX_INDEX + 1), repeat=3):
            if not any(index) or _primitive(index) != index:
                continue
            trial = self.cartesian(np.array(index, float))
            cos = abs(float(unit @ trial)) / np.linalg.norm(trial)
            angle = float(np.degrees(np.arccos(min(cos, 1.0))))
            # Prefer the simpler index unless the other is clearly closer
            if angle < best_angle - 0.5 or (
                angle < best_angle + 0.5 and sum(map(abs, index)) < sum(map(abs, best))
            ):
                best, best_angle = index, angle
        axes = self.direct / np.linalg.norm(self.direct, axis=1, keepdims=True)
        axis_angles = tuple(
            float(np.degrees(np.arccos(min(abs(float(unit @ axis)), 1.0)))) for axis in axes
        )
        return Direction(np.asarray(vector, float), best, best_angle, axis_angles)

    def dimer_vector(self, coords: np.ndarray, molecule: np.ndarray) -> np.ndarray | None:
        """Fractional centroid-to-centroid vector of a dimer file, or None if
        the file's atoms are not in the CIF's Cartesian frame."""
        ends = []
        for m in (0, 1):
            centre = self.fractional(coords[molecule == m].mean(axis=0))
            gaps = [
                self._length((centre - site) - np.round(centre - site)) for _, site in self._sites
            ]
            if min(gaps) > FRAME_TOLERANCE:
                return None
            ends.append(centre)
        return ends[1] - ends[0]

    # -- motifs ---------------------------------------------------------------

    def _site(self, label: str, point: np.ndarray) -> tuple[tuple, np.ndarray]:
        """(site key, lattice translation) of a centroid."""
        rounded = np.round(point, 4)
        shift = np.floor(rounded)
        base = np.round(rounded - shift, 3) % 1.0
        return (label, *base.tolist()), shift.astype(int)

    def motif(self, pairs: list[Pair]) -> Motif:
        """What the interaction(s) in ``pairs`` build when copied by every symmetry operation."""
        edges: set[tuple] = set()
        for pair in pairs:
            ends = (pair.first_op(self.centroids[pair.first]),
                    pair.second_op(self.centroids[pair.second]))
            for s in self.symops:
                (k1, t1), (k2, t2) = (self._site(pair.first, s(ends[0])),
                                      self._site(pair.second, s(ends[1])))
                shift = tuple((t2 - t1).tolist())
                edges.add((k1, k2, shift))
                edges.add((k2, k1, tuple(-v for v in shift)))

        neighbours: dict[tuple, list[tuple]] = {}
        for k1, k2, shift in edges:
            neighbours.setdefault(k1, []).append((k2, np.array(shift)))

        start = self._site(pairs[0].first, pairs[0].first_op(self.centroids[pairs[0].first]))[0]
        offset = {start: np.zeros(3, int)}
        queue = deque([start])
        cycles = []
        while queue:
            node = queue.popleft()
            for other, shift in neighbours.get(node, []):
                reached = offset[node] + shift
                if other not in offset:
                    offset[other] = reached
                    queue.append(other)
                elif np.any(reached != offset[other]):
                    cycles.append(reached - offset[other])

        rank = int(np.linalg.matrix_rank(np.array(cycles, float))) if cycles else 0
        if rank == 0:
            # Each copy of a finite motif is the component itself
            degrees = [len({(o, tuple(s)) for o, s in neighbours.get(n, [])}) for n in offset]
            n = len(offset)
            return Motif(0, n_molecules=n, ring=n > 2 and all(d == 2 for d in degrees))
        if rank == 1:
            direction = np.array(_primitive(cycles[0]))
            steps = [abs(round(float(c @ direction) / float(direction @ direction))) for c in cycles]
            generator = reduce(gcd, steps) * direction
            return Motif(1, index=tuple(direction.tolist()),
                         repeat=self._length(generator.astype(float)))
        if rank == 2:
            first = cycles[0]
            second = next(c for c in cycles[1:] if np.any(np.cross(first, c)))
            return Motif(2, index=_primitive(np.cross(first, second)))
        return Motif(3)


def find_cif(dimers_folder: Path, input_dir: Path | None = None) -> Path | None:
    """The CIF of the crystal: ``<name>.cif`` beside ``<name>_dimers``, else a
    lone CIF there, else the shallowest CIF under ``input_dir``."""
    dimers_folder = Path(dimers_folder)
    stem = re.sub(r"[_\-.]?dimers$", "", dimers_folder.name, flags=re.IGNORECASE)
    parent = dimers_folder.parent
    if stem and (parent / f"{stem}.cif").exists():
        return parent / f"{stem}.cif"
    nearby = sorted(parent.glob("*.cif"))
    if len(nearby) == 1:
        return nearby[0]
    matching = [p for p in nearby if stem and p.stem.lower().startswith(stem.lower())]
    if matching:
        return matching[0]
    if input_dir is not None and Path(input_dir).is_dir():
        found = sorted(Path(input_dir).rglob("*.cif"), key=lambda p: (len(p.parts), str(p)))
        if found:
            return found[0]
    return nearby[0] if nearby else None


__all__ = ["CrystalMotifs", "Direction", "Motif", "Op", "Pair", "find_cif"]
