from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Callable, Iterable, Union

import numpy as np
import trimesh

LOG = logging.getLogger("SHAPECLASS")


@dataclass
class Frame:
    """Single crystal shape frame. Stores raw particle data and exposes named column slices.

    Column layout (CrystalGrower XYZ format):
        0  - particle type
        1  - particle number
        2  - layer  (values >= 99 are sentinel values)
        3  - x
        4  - y
        5  - z
        6  - site number  (newer CG versions only)
        7  - particle energy  (newer CG versions only)
    """

    raw: np.ndarray
    comment: Optional[str] = None
    validate_cg: bool = False

    def __post_init__(self) -> None:
        if self.validate_cg and (self.raw.ndim != 2 or self.raw.shape[1] < 6):
            LOG.warning(
                "Frame does not appear to be a valid CrystalGrower format: "
                "expected at least 6 columns, got shape %s. "
                "Column-based properties (layer, particle_type, etc.) may not be available.",
                self.raw.shape,
            )

    # --- named column properties ---
    @property
    def particle_type(self) -> np.ndarray:
        return self.raw[:, 0]

    @property
    def number(self) -> np.ndarray:
        return self.raw[:, 1]

    @property
    def layer(self) -> np.ndarray:
        return self.raw[:, 2]

    @property
    def coords(self) -> np.ndarray:
        """Spatial coordinates (x, y, z). Falls back to first 3 columns for non-CG formats."""
        return self.raw[:, 3:6] if self.raw.shape[1] >= 6 else self.raw[:, :3]

    @property
    def x(self) -> np.ndarray:
        return self.coords[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.coords[:, 1]

    @property
    def z(self) -> np.ndarray:
        return self.coords[:, 2]

    @property
    def site(self) -> Optional[np.ndarray]:
        return self.raw[:, 6] if self.raw.shape[1] >= 7 else None

    @property
    def energy(self) -> Optional[np.ndarray]:
        return self.raw[:, 7] if self.raw.shape[1] >= 8 else None

    # --- container behaviour ---
    def __len__(self) -> int:
        return len(self.raw)

    def __getitem__(self, idx: Union[int, slice]) -> np.ndarray:
        return self.coords[idx]

    def __iter__(self) -> Iterable[np.ndarray]:
        return iter(self.coords)


@dataclass
class Frames:
    """Container for multiple frames. Behaves like a list of Frame objects."""

    _frames: list[Frame] = field(default_factory=list)

    # --- core list-like behaviour ---
    def __len__(self) -> int:
        return len(self._frames)

    def __getitem__(self, idx: Union[int, slice]) -> Union[Frame, "Frames"]:
        if isinstance(idx, slice):
            return Frames(self._frames[idx])
        return self._frames[idx]

    def __iter__(self) -> Iterable[Frame]:
        return iter(self._frames)

    def append(self, frame: Frame) -> None:
        self._frames.append(frame)

    def extend(self, frames: Iterable[Frame]) -> None:
        self._frames.extend(frames)

    # --- convenience views ---
    @property
    def coords(self) -> dict[int, np.ndarray]:
        """All frame coordinates as dict {index: coords}."""
        return {i: f.coords for i, f in enumerate(self._frames)}

    @property
    def comments(self) -> dict[int, Optional[str]]:
        """All frame comments as dict {index: comment}."""
        return {i: f.comment for i, f in enumerate(self._frames)}

    def get_coords(self, idx: int) -> Optional[np.ndarray]:
        """Convenience: coords for a single frame."""
        if 0 <= idx < len(self._frames):
            return self._frames[idx].coords
        return None


# ---------------- CrystalShape ----------------
@dataclass
class CrystalShape:
    """Handles crystal shape data from various file formats."""

    filepath: Path
    frames: Frames = field(default_factory=Frames)
    xyz: Optional[np.ndarray] = None

    # ---- Core container behaviour ----
    def __len__(self) -> int:
        return len(self.frames)

    def __getitem__(self, idx: int) -> Frame:
        return self.frames[idx]

    def __iter__(self):
        return iter(self.frames)

    # ---- Convenience views ----
    @property
    def movie(self) -> dict[int, np.ndarray]:
        """Return all frames as dict {index: coords}."""
        return self.frames.coords

    @property
    def coords(self) -> Optional[np.ndarray]:
        """Return the coordinates of the last frame (index -1)."""
        return self.frames.get_coords(-1)

    # ---- Parsing helpers ----
    @staticmethod
    def parse_xyz_file(
        filepath: Path,
        progress_callback: Optional[Callable[[int, int], None]] = None,
        clean: bool = True,
    ) -> Frames:
        """Parse multi-frame XYZ into Frames container."""
        frames = Frames()

        with filepath.open("r", encoding="utf-8") as file:
            frame_idx = 0
            while True:
                header = file.readline()
                if not header:
                    break

                try:
                    n_atoms = int(header.strip())
                except ValueError as e:
                    raise ValueError(f"Invalid XYZ header at frame {frame_idx}: {e}")

                comment = file.readline().strip()

                try:
                    raw = np.loadtxt(file, max_rows=n_atoms, dtype=float, ndmin=2)
                except ValueError as e:
                    if clean:
                        raw_text = Path(filepath).read_text(encoding="utf-8").replace("*", "0")
                        Path(filepath).write_text(raw_text, encoding="utf-8")
                        return CrystalShape.parse_xyz_file(filepath, progress_callback, clean=False)
                    raise e

                frames.append(Frame(raw=raw, comment=comment, validate_cg=True))
                frame_idx += 1

                if progress_callback:
                    try:
                        total_frames = int(comment.split("//")[1])
                    except Exception:
                        total_frames = frame_idx
                    progress_callback(frame_idx, total_frames)

        return frames

    @staticmethod
    def normalise_verts(verts, center=True):
        if center:
            verts = verts - np.mean(verts, axis=0)
        norm = np.linalg.norm(verts, axis=1).max()
        verts /= norm
        return verts

    # ---- Loader ----
    @classmethod
    def from_file(
        cls,
        filepath: Path,
        progress_callback: Optional[Callable[[int, int], None]] = None,
        normalise: bool = True,
    ) -> "CrystalShape":
        """Factory method to create CrystalShape from .xyz, .txt, .stl, .glb, .ply."""
        filepath = Path(filepath)
        suffix = filepath.suffix.lower()

        if suffix == ".xyz":
            frames = cls.parse_xyz_file(filepath, progress_callback)
            xyz = frames.get_coords(0)

        elif suffix == ".txt":
            arr = np.loadtxt(filepath, skiprows=2)
            frames = Frames([Frame(raw=arr, comment="txt-file")])
            xyz = frames[0].coords

        elif suffix in {".stl", ".glb", ".ply"}:
            mesh = trimesh.load(filepath)
            frames = Frames([Frame(raw=mesh.vertices, comment=suffix)])
            xyz = frames[0].coords

        else:
            raise ValueError(f"Unsupported file format: {suffix}")

        if xyz is not None and normalise:
            xyz = cls.normalise_verts(xyz)

        return cls(filepath=filepath, frames=frames, xyz=xyz)

    def get_frame_coords(self, frame_idx: int = 0) -> Optional[np.ndarray]:
        return self.frames.get_coords(frame_idx)

    def get_all_frame_coords(self) -> dict[int, np.ndarray]:
        return self.frames.coords

    @staticmethod
    def write_xyz(frames: Frames, out_path: Path) -> None:
        """Write a Frames container to a multi-frame XYZ file."""
        with out_path.open("w", encoding="utf-8") as f:
            for frame in frames:
                f.write(f"{len(frame)}\n")
                f.write(f"{frame.comment or ''}\n")
                np.savetxt(f, frame.raw, fmt="%.6g")

    @staticmethod
    def count_frames(filepath: Path) -> int:
        """Count frames by scanning headers only — does not load atom data."""
        count = 0
        with filepath.open("r", encoding="utf-8") as f:
            while True:
                header = f.readline()
                if not header:
                    break
                try:
                    n_atoms = int(header.strip())
                except ValueError:
                    break
                f.readline()  # comment line
                for _ in range(n_atoms):
                    f.readline()
                count += 1
        return count

    @staticmethod
    def iter_frames(filepath: Path):
        """Yield (index, Frame) one at a time — constant memory regardless of file size."""
        with filepath.open("r", encoding="utf-8") as f:
            idx = 0
            while True:
                header = f.readline()
                if not header:
                    break
                try:
                    n_atoms = int(header.strip())
                except ValueError:
                    break
                comment = f.readline().strip()
                raw = np.loadtxt(f, max_rows=n_atoms, dtype=float, ndmin=2)
                yield idx, Frame(raw=raw, comment=comment)
                idx += 1


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        prog="crystal_io",
        description="Inspect and manipulate CrystalGrower XYZ files.",
    )
    parser.add_argument("input", type=Path, help="Input XYZ file")

    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("info", help="Print frame count and atoms-per-frame summary")

    p_ex = sub.add_parser(
        "extract",
        help="Extract frame(s): 1 arg=single frame (negative ok); 2=start stop; 3=start stop step",
    )
    p_ex.add_argument("nums", nargs="+", type=int, metavar="N")
    p_ex.add_argument("-o", "--output", type=Path, required=True, help="Output XYZ file")

    args = parser.parse_args()

    if args.cmd == "info":
        n = 0
        min_atoms = max_atoms = total_atoms = 0
        for i, frame in CrystalShape.iter_frames(args.input):
            size = len(frame)
            if i == 0:
                min_atoms = max_atoms = size
            else:
                min_atoms = min(min_atoms, size)
                max_atoms = max(max_atoms, size)
            total_atoms += size
            n += 1
        print(f"Frames : {n}")
        if n:
            print(f"Atoms  : min={min_atoms}  max={max_atoms}  mean={total_atoms/n:.1f}")

    elif args.cmd == "extract":
        nums = args.nums
        if len(nums) not in (1, 2, 3):
            print("Error: extract takes 1, 2, or 3 integers", file=sys.stderr)
            sys.exit(1)

        # resolve negative indices without loading frame data
        if any(x < 0 for x in nums):
            n = CrystalShape.count_frames(args.input)
            nums = [x if x >= 0 else n + x for x in nums]

        if len(nums) == 1:
            target = {nums[0]}
        else:
            start, stop = nums[0], nums[1]
            step = nums[2] if len(nums) == 3 else 1
            target = set(range(start, stop, step))

        if not target:
            print("No frames in selection.", file=sys.stderr)
            sys.exit(1)

        max_target = max(target)
        count = 0
        with args.output.open("w", encoding="utf-8") as out_f:
            for idx, frame in CrystalShape.iter_frames(args.input):
                if idx > max_target:
                    break
                if idx in target:
                    out_f.write(f"{len(frame)}\n")
                    out_f.write(f"{frame.comment or ''}\n")
                    np.savetxt(out_f, frame.raw, fmt="%.6g")
                    count += 1

        print(f"Wrote {count} frame(s) → {args.output}")
