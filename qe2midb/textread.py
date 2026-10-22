"""Readers of the text outputs of QE post-processing programs.

dos.x (fildos), projwfc.x (filpdos files and stdout), matdyn.x (flfrq, flfrq.gp, fldos),
dynmat.x (filout and stdout), ph.x (dynamical matrix files, input), pw.x input (K_POINTS
labels). Tolerant to spin-polarized, non-collinear and spin-orbit layouts and to numbered
species labels (Fe1, Fe2).
"""
import glob
import os
import re

import numpy as np

from qe2midb.xmlread import element_of

FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"


def head(path, nbytes=4096):
    try:
        with open(path, "rb") as f:
            return f.read(nbytes).decode("utf-8", "replace")
    except OSError:
        return ""


def program_version(path):
    """(PROGRAM, version) from the 'Program XXX v.Y starts' banner of a QE stdout."""
    m = re.search(r"Program (\S+) v\.(\S+) starts", head(path, 20000))
    return (m.group(1), m.group(2)) if m else (None, None)


def _loadtxt(path):
    return np.loadtxt(path, comments="#", ndmin=2)


# ---------------------------------------------------------------------- dos.x
def is_dos(path):
    h = head(path, 300).lstrip()
    return h.startswith("#") and "E (eV)" in h and "dos" in h.lower() and "EFermi" in h


def dos_file(path):
    """dos.x output: energy (eV), DOS [spin][E], integrated DOS, E_F from the header."""
    first = head(path, 300).splitlines()[0]
    m = re.search(r"EFermi\s*=\s*(" + FLOAT + ")", first)
    data = _loadtxt(path)
    dos = data[:, 1:3].T if "dosup(E)" in first.replace("#", " ").split() else data[:, 1:2].T
    return {"energy": data[:, 0], "total": dos, "integrated": data[:, -1],
            "efermi": float(m.group(1)) if m else None}


# ---------------------------------------------------------------------- projwfc.x
_PDOS_RE = re.compile(r"\.pdos_atm#(\d+)\(([^)]+)\)_wfc#(\d+)\(([^)]+)\)$")


def pdos_prefixes(directory):
    """Prefixes (path without '.pdos_tot') of the projwfc.x PDOS sets in a directory."""
    return sorted(p[: -len(".pdos_tot")] for p in glob.glob(os.path.join(directory, "*.pdos_tot")))


def _ldos_columns(header, data):
    """Per-wavefunction file -> local DOS [channel][E]."""
    if header.lstrip("# ").startswith("ik"):
        raise ValueError("k-resolved PDOS (kresolveddos=.true.) is not supported")
    if "ldosup(E)" in header.replace("#", " ").split():   # LSDA or non-collinear without spin-orbit
        return data[:, 1:3].T
    return data[:, 1:2].T


def pdos_files(prefix, species_map=None, noncolin=False):
    """PDOS of projwfc.x summed per (atom, orbital l).

    prefix: path without the '.pdos_*' suffix (e.g. 'out/pw'). Orbitals with the same l on
    the same atom (semicore + valence, or j = l±1/2 with spin-orbit) are summed. Species
    labels are mapped to elements with species_map (from the XML) when given. In
    non-collinear runs the two spin projections are summed (one channel)."""
    energy, acc = None, {}
    for f in sorted(glob.glob(glob.escape(prefix) + ".pdos_atm#*")):
        m = _PDOS_RE.search(f)
        if not m:
            continue
        iat, label, orb = int(m.group(1)), m.group(2), m.group(4)
        l = orb.split("_")[0][0]                 # 's', 'p_j1.5' -> 'p', '5d' unlikely
        with open(f) as fh:
            header = fh.readline()
        data = _loadtxt(f)
        ld = _ldos_columns(header, data)
        if noncolin and ld.shape[0] == 2:
            ld = ld.sum(axis=0, keepdims=True)
        energy = data[:, 0]
        key = (iat, label, l)
        acc[key] = acc[key] + ld if key in acc else ld
    if not acc:
        raise FileNotFoundError(f"no PDOS files for prefix {prefix}")
    smap = species_map or {}
    entries = [{"atom_index": k[0], "atom": f"{smap.get(k[1], element_of(k[1]))}{k[0]}",
                "species": smap.get(k[1], element_of(k[1])), "label": k[1], "orbital": k[2], "values": v}
               for k, v in sorted(acc.items())]
    total = pdos_sum = None
    tot = prefix + ".pdos_tot"
    if os.path.exists(tot):
        with open(tot) as fh:
            header = fh.readline()
        data = _loadtxt(tot)
        cols = header.replace("#", " ").split()
        if "dosup(E)" in cols:                   # LSDA: dosup dosdw pdosup pdosdw
            total, pdos_sum = data[:, 1:3].T, data[:, 3:5].T
        elif "pdosup(E)" in cols:                # non-collinear: dos, pdosup, pdosdw
            total, pdos_sum = data[:, 1:2].T, data[:, 2:4].sum(axis=1)[None]
        else:
            total, pdos_sum = data[:, 1:2].T, data[:, 2:3].T
    return {"energy": energy, "entries": entries, "total": total, "pdos_sum": pdos_sum}


def lowdin_spilling(projwfc_stdout):
    m = re.search(r"Spilling Parameter:\s*(" + FLOAT + ")", open(projwfc_stdout, errors="replace").read())
    return float(m.group(1)) if m else None


# ---------------------------------------------------------------------- namelists / inputs
def namelist(text, name):
    """Very small Fortran namelist reader: {key(lowercase): raw value string}."""
    m = re.search(r"&" + name + r"\b(.*?)\n\s*/", text, re.S | re.I)
    if not m:
        return None
    body = re.sub(r"!.*", "", m.group(1))
    out = {}
    for k, v in re.findall(r"([A-Za-z_][\w()%]*)\s*=\s*('[^']*'|\"[^\"]*\"|[^,\n]+)", body):
        out[k.lower()] = v.strip().strip("'\"")
    return out


def fbool(v):
    return str(v).strip().lower() in (".true.", "true", "t", ".t.")


def qe_input_kind(path):
    """'pw', 'ph', 'matdyn', 'dynmat', 'q2r', 'dos', 'projwfc', 'bands' or None (by namelist)."""
    t = head(path, 20000).lower()
    if "&control" in t and ("&system" in t or "atomic_positions" in t):
        return "pw"
    if "&inputph" in t:
        return "ph"
    if "&projwfc" in t:
        return "projwfc"
    if "&dos" in t:
        return "dos"
    if "&bands" in t:
        return "bands"
    if "&input" in t:
        has_fildyn, has_flfrc = "fildyn" in t, "flfrc" in t
        if "flfrq" in t or "fldos" in t or (has_flfrc and not has_fildyn):
            return "matdyn"
        if has_fildyn and has_flfrc:
            return "q2r"
        if has_fildyn:
            return "dynmat"
    return None


def pw_kpoint_labels(pw_input_path):
    """Labels of a K_POINTS crystal_b / tpiba_b card: [(frac or cart coords, label)] or []."""
    text = open(pw_input_path, errors="replace").read()
    m = re.search(r"K_POINTS\s*[({]?\s*(crystal_b|tpiba_b)\s*[)}]?\s*\n\s*(\d+)\s*\n", text, re.I)
    if not m:
        return []
    n = int(m.group(2))
    lines = text[m.end():].splitlines()[:n]
    out = []
    for line in lines:
        nums = re.findall(FLOAT, line.split("!")[0])
        label = line.split("!")[1].strip() if "!" in line else ""
        if len(nums) >= 3:
            out.append(([float(x) for x in nums[:3]], label))
    return out


def pw_input(path):
    """Main parameters of a pw.x input: calculation, occupations/smearing/degauss (Ry),
    K_POINTS automatic grid."""
    text = open(path, errors="replace").read()
    ctl = namelist(text, "control") or {}
    sysn = namelist(text, "system") or {}
    out = {"calculation": ctl.get("calculation", "scf").lower(),
           "occupations": sysn.get("occupations", "fixed").lower()}
    if out["occupations"] == "smearing":
        out["smearing"] = {"occupations": "smearing", "smearing": sysn.get("smearing", "gaussian"),
                           "degauss": float(sysn.get("degauss", "0").lower().replace("d", "e"))}
    m = re.search(r"K_POINTS\s*[({]?\s*automatic\s*[)}]?\s*\n\s*(\d+)\s+(\d+)\s+(\d+)", text, re.I)
    out["kgrid"] = [int(m.group(i)) for i in (1, 2, 3)] if m else None
    if re.search(r"K_POINTS\s*[({]?\s*gamma", text, re.I):
        out["kgrid"] = [1, 1, 1]
    return out


def pw_stdout_scf_energy(path):
    """Final total energy (eV) of a pw.x stdout of a converged single-point scf run, else None
    (relax/md runs, which print 'Begin final coordinates' or several ionic steps, are refused)."""
    t = open(path, errors="replace").read()
    if "End of self-consistent calculation" not in t or "Begin final coordinates" in t \
            or "BFGS Geometry Optimization" in t or "Molecular Dynamics Calculation" in t:
        return None
    if t.count("Program PWSCF") != 1 or t.count("convergence has been achieved") != 1:
        return None                              # several runs appended, or not converged
    e = re.findall(r"^!\s+total energy\s+=\s+(" + FLOAT + r")\s+Ry", t, re.M)
    return float(e[-1]) * 13.605693122994 if e else None   # DFT+U may print '!' every iteration


def ph_input(path):
    """ph.x input: qgrid (ldisp) or Gamma only, epsil, fildyn."""
    nl = namelist(open(path, errors="replace").read(), "inputph") or {}
    q = [int(nl.get(f"nq{i}", 1)) for i in (1, 2, 3)] if fbool(nl.get("ldisp", "false")) else [1, 1, 1]
    return {"qgrid": q, "ldisp": fbool(nl.get("ldisp", "false")), "epsil": fbool(nl.get("epsil", "false")),
            "fildyn": nl.get("fildyn", "matdyn")}


# ---------------------------------------------------------------------- ph.x dynamical matrices
def is_dyn(path):
    return "Dynamical matrix file" in head(path, 200)


def dyn0_grid(path):
    """q-grid from the dyn0 file written by ph.x with ldisp (first line: nq1 nq2 nq3)."""
    first = head(path, 200).splitlines()
    try:
        q = [int(x) for x in first[0].split()[:3]]
        return q if len(q) == 3 else None
    except (ValueError, IndexError):
        return None


def epsilon_zstar(dyn_path):
    """High-frequency dielectric tensor and Born charges Z*[atom][3][3] (E-U block) of a
    Gamma dynamical-matrix file; (None, []) when absent."""
    txt = open(dyn_path, errors="replace").read()
    eps, zs = None, []
    m = re.search(r"Dielectric Tensor:\s*\n((?:\s*" + FLOAT + r"\s+" + FLOAT + r"\s+" + FLOAT + r"\s*\n){3})", txt)
    if m:
        eps = [[float(x) for x in l.split()] for l in m.group(1).strip().splitlines()]
    m = re.search(r"Effective Charges E-U.*?\n(.*?)(?=Effective Charges U-E|Diagonalizing|\Z)", txt, re.S)
    if m:
        for blk in re.split(r"atom #\s*\d+", m.group(1))[1:]:
            vals = [float(x) for x in re.findall(r"[-]?\d+\.\d+", blk)][:9]
            if len(vals) == 9:
                zs.append(np.array(vals).reshape(3, 3).tolist())
    return eps, zs


def is_gamma_dyn(path):
    m = re.search(r"Diagonalizing the dynamical matrix\s*\n\s*q = \(\s*(" + FLOAT + r")\s+(" + FLOAT + r")\s+("
                  + FLOAT + r")\s*\)", open(path, errors="replace").read())
    return bool(m) and all(abs(float(m.group(i))) < 1e-8 for i in (1, 2, 3))


# ---------------------------------------------------------------------- modes (dynmat.x / matdyn.x)
_FREQ_RE = re.compile(r"freq\s*\(\s*\d+\)\s*=\s*(" + FLOAT + r")\s*\[THz\]\s*=\s*(" + FLOAT + r")\s*\[cm-1\]")


def modes_file(path):
    """Normal modes file (dynmat.x filout or matdyn.x flvec): list of q blocks, each
    {q, modes: [{freq_THz, freq_cm1, eigvec [nat][3] (real part)}]}."""
    lines = open(path, errors="replace").read().splitlines()
    blocks, cur, i = [], None, 0
    while i < len(lines):
        line = lines[i]
        mq = re.match(r"\s*q\s*=\s*(" + FLOAT + r")\s+(" + FLOAT + r")\s+(" + FLOAT + r")", line)
        if mq:
            cur = {"q": [float(mq.group(k)) for k in (1, 2, 3)], "modes": []}
            blocks.append(cur)
        m = _FREQ_RE.search(line)
        if m and cur is not None:
            vec = []
            i += 1
            while i < len(lines) and lines[i].strip().startswith("("):
                nums = [float(x) for x in re.findall(FLOAT, lines[i])]
                vec.append(nums[0::2][:3])
                i += 1
            cur["modes"].append({"freq_THz": float(m.group(1)), "freq_cm1": float(m.group(2)), "eigvec": vec})
            continue
        i += 1
    return blocks


def is_modes_file(path):
    h = head(path, 3000)
    return "diagonalizing the dynamical matrix" in h.lower() and "freq (" in h


def dynmat_ir(stdout_path):
    """Final table of dynmat.x stdout: mode, frequency (cm-1), IR intensity ((D/A)^2/amu)."""
    txt = open(stdout_path, errors="replace").read()
    if "# mode" not in txt:
        return []
    rows = []
    for line in txt.split("# mode")[-1].splitlines()[1:]:
        p = line.split()
        if len(p) >= 3 and p[0].isdigit():
            rows.append({"mode": int(p[0]), "freq_cm1": float(p[1]),
                         "ir": float(p[3]) if len(p) >= 4 else None})
    return rows


# ---------------------------------------------------------------------- matdyn.x
def matdyn_input(path):
    """Namelist + q-point list of a matdyn.x input."""
    text = open(path, errors="replace").read()
    nl = namelist(text, "input") or {}
    qpts = []
    tail = text[re.search(r"&input.*?\n\s*/", text, re.S | re.I).end():] if re.search(r"&input.*?\n\s*/", text, re.S | re.I) else ""
    lines = [l for l in tail.splitlines() if l.strip()]
    if lines and not fbool(nl.get("dos", "false")):
        try:
            n = int(lines[0].split()[0])
            for line in lines[1:1 + n]:
                nums = re.findall(FLOAT, line.split("!")[0])
                label = line.split("!")[1].strip() if "!" in line else ""
                qpts.append(([float(x) for x in nums[:3]], int(float(nums[3])) if len(nums) > 3 else 1, label))
        except (ValueError, IndexError):
            qpts = []
    return {"dos": fbool(nl.get("dos", "false")), "flfrq": nl.get("flfrq", "matdyn.freq"),
            "fldos": nl.get("fldos", "matdyn.dos"), "flvec": nl.get("flvec", "matdyn.modes"),
            "q_in_band_form": fbool(nl.get("q_in_band_form", "false")),
            "q_in_cryst_coord": fbool(nl.get("q_in_cryst_coord", "false")), "qpoints": qpts}


def matdyn_freq(path):
    """matdyn flfrq (plot format): q-points (cartesian 2π/alat) and frequencies [q][mode] (cm-1).

    Each block is one line with the q-point (3 coordinates; QE >= 7 may append a 4th
    number) followed by lines with the nbnd frequencies."""
    lines = open(path, errors="replace").read().splitlines()
    m = re.search(r"&plot\s+nbnd\s*=\s*(\d+)\s*,\s*nks\s*=\s*(\d+)", lines[0] if lines else "")
    if not m:
        raise ValueError(f"{path}: not a matdyn frequency file")
    nb, nk = int(m.group(1)), int(m.group(2))
    q, w, i = [], [], 1
    while len(q) < nk and i < len(lines):
        nums = [float(x) for x in re.findall(FLOAT, lines[i])]
        i += 1
        if len(nums) < 3:
            continue
        q.append(nums[:3])
        vals = []
        while len(vals) < nb and i < len(lines):
            vals += [float(x) for x in re.findall(FLOAT, lines[i])]
            i += 1
        w.append(vals[:nb])
    return np.array(q), np.array(w)


def plot_file_nbnd(path):
    m = re.search(r"&plot\s+nbnd\s*=\s*(\d+)", head(path, 200))
    return int(m.group(1)) if m else None


def is_matdyn_freq(path):
    return bool(re.search(r"&plot\s+nbnd\s*=", head(path, 200)))


def is_phdos(path):
    return "Frequency[cm^-1]" in head(path, 200).replace(" ", "") or "Frequency[cm-1]" in head(path, 200).replace(" ", "")


def phdos(path, nat=None):
    """matdyn fldos: frequency (cm-1), total DOS, per-atom projections when present."""
    data = _loadtxt(path)
    per_atom = data[:, 2:].T if data.shape[1] > 2 else None
    if per_atom is not None and nat is not None and per_atom.shape[0] != nat:
        per_atom = None                          # older matdyn: single PDOS column
    return {"freq": data[:, 0], "total": data[:, 1], "per_atom": per_atom}
