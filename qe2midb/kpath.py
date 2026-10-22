"""Band-structure path from the k-points themselves (no ASE BandPath file needed).

Given k-points in crystal coordinates and the cell, builds the linear axis (1/Å, 2π
included), detects discontinuities (a step much larger than the typical step; the
convention of the normalized data is a repeated x; local rule against both neighbouring
steps), finds the path vertices (direction
changes) and labels them:
  1. labels given explicitly (K_POINTS crystal_b comments, matdyn input), else
  2. high-symmetry points of the lattice from ASE (if installed), else pymatgen, else
  3. no label ('').
"""
import numpy as np

GAMMA_NAMES = {"G", "GAMMA", "\\GAMMA", "Γ", "GM"}
JUMP_FACTOR = 4.0      # a step longer than JUMP_FACTOR x both neighbouring steps is a discontinuity
ANGLE_TOL = 1e-3       # radians; direction change that defines a vertex


def reciprocal(cell):
    """Reciprocal vectors (rows, 1/Å, with 2π) of a cell given as rows in Å."""
    return 2 * np.pi * np.linalg.inv(np.asarray(cell, dtype=float)).T


def pretty(label):
    if label is None:
        return ""
    return "Γ" if label.strip().upper() in GAMMA_NAMES else label.strip()


def special_points(cell):
    """{label: fractional coords} of the high-symmetry points of the lattice, or {}."""
    try:
        from ase.cell import Cell
        return {k: np.asarray(v, dtype=float) for k, v in Cell(np.asarray(cell)).bandpath(npoints=0).special_points.items()}
    except Exception:  # noqa: BLE001 - optional dependency or unusual cell
        pass
    try:
        from pymatgen.core import Lattice, Structure
        from pymatgen.symmetry.bandstructure import HighSymmKpath
        st = Structure(Lattice(cell), ["H"], [[0, 0, 0]])
        kp = HighSymmKpath(st, path_type="setyawan_curtarolo")
        if np.allclose(kp.prim.lattice.matrix, st.lattice.matrix, atol=1e-6):
            return {k: np.asarray(v, dtype=float) for k, v in kp.kpath["kpoints"].items()}
    except Exception:  # noqa: BLE001
        pass
    return {}


def _match(k, points, tol=1e-4):
    for label, p in points.items():
        d = np.asarray(k) - p
        if np.all(np.abs(d - np.round(d)) < tol):
            return label
    return None


def build(kfrac, cell, labels=None):
    """Linear axis, ticks and path string of a band/phonon path.

    kfrac: [n][3] crystal coordinates; labels: optional {label: frac coords} with priority
    over the lattice special points. Returns {x, ticks, path, kpoints_frac, n_labeled}."""
    kfrac = np.asarray(kfrac, dtype=float)
    kc = kfrac @ reciprocal(cell)
    steps = np.linalg.norm(np.diff(kc, axis=0), axis=1)
    # a discontinuity is a step much longer than BOTH neighbouring steps (local rule: paths
    # with a fixed number of points per segment have different step lengths per segment)
    jump = np.zeros_like(steps, dtype=bool)
    for i, st in enumerate(steps):
        neigh = [steps[j] for j in (i - 1, i + 1) if 0 <= j < len(steps) and steps[j] > 1e-9]
        jump[i] = bool(neigh) and st > JUMP_FACTOR * max(neigh)
    x = np.concatenate([[0.0], np.cumsum(np.where(jump, 0.0, steps))])

    vertices = {0, len(kfrac) - 1}
    for j in np.where(jump)[0]:
        vertices.update((int(j), int(j) + 1))
    for i in range(1, len(kfrac) - 1):
        if jump[i - 1] or jump[i]:
            continue
        v1, v2 = kc[i] - kc[i - 1], kc[i + 1] - kc[i]
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-9 or n2 < 1e-9:
            continue
        if np.arccos(np.clip(np.dot(v1, v2) / (n1 * n2), -1, 1)) > ANGLE_TOL:
            vertices.add(i)

    points = {}
    if labels:
        points.update({k: np.asarray(v, dtype=float) for k, v in labels.items() if k})
    lattice_points = special_points(cell)
    ticks, path, n_labeled = [], "", 0
    for i in sorted(vertices):
        label = _match(kfrac[i], points) or _match(kfrac[i], lattice_points)
        n_labeled += label is not None
        ticks.append({"x": round(float(x[i]), 6), "label": pretty(label)})
        if i > 0 and jump[i - 1]:
            path += ","
        path += "G" if pretty(label) == "Γ" else (label or "?")
    return {"x": x, "ticks": ticks, "path": path, "kpoints_frac": kfrac, "n_labeled": n_labeled,
            "n_vertices": len(vertices)}
