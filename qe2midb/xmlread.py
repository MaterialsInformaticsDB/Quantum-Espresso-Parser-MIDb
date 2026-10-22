"""Reader of the pw.x XML data file (data-file-schema.xml / <prefix>.xml), any QE >= 6.3.

The schema declared by the file (xsi:schemaLocation) is looked up among the official
schemas vendored in ``schemas/`` (QEF/qeschemas); when absent, the schema with the
nearest date is used. Decoding is lenient ("lax"): QE itself writes files with small
deviations from its own schema, and the values do not depend on which recent schema is
used (checked on QE 6.3-7.6 against the text output). The schema used and the number of
validation errors are kept for provenance.
"""
import glob
import os
import re
import warnings

import numpy as np
import xmlschema

HARTREE_EV = 27.211386245988
BOHR_ANG = 0.529177210903
SCHEMA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "schemas")
MIN_VERSION = (6, 3)

_SCHEMAS = {}


class QEXMLError(Exception):
    pass


def _schema(path):
    if path not in _SCHEMAS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _SCHEMAS[path] = xmlschema.XMLSchema(path)
    return _SCHEMAS[path]


def _schema_date(name):
    """qes_YYMMDD / qes_dev_YYMMDD / qes-YYYYMMDD -> sortable int (YYYYMMDD)."""
    m = re.search(r"(\d{6,8})", name or "")
    if not m:
        return None
    s = m.group(1)
    return int(s) if len(s) == 8 else 20000000 + int(s)


def available_schemas():
    return sorted(glob.glob(os.path.join(SCHEMA_DIR, "qes*.xsd")))


def schema_for(declared):
    """Path of the vendored schema for a declared name (exact, dev variant or nearest date)."""
    if declared:
        for cand in (declared, declared.replace("qes_", "qes_dev_"), declared.replace("-", "_")):
            p = os.path.join(SCHEMA_DIR, cand + ".xsd")
            if os.path.exists(p):
                return p
    target = _schema_date(declared)
    dated = [(d, p) for p in available_schemas() if (d := _schema_date(os.path.basename(p))) is not None]
    if not dated:
        raise QEXMLError("no vendored QE schema found")
    if target is None:
        return max(dated)[1]
    return min(dated, key=lambda dp: abs(dp[0] - target))[1]


def sniff(path, nbytes=6000):
    """(is_pw_xml, declared_schema, qe_version) from the head of a file, without parsing."""
    try:
        with open(path, "rb") as f:
            head = f.read(nbytes).decode("utf-8", "replace")
    except OSError:
        return False, None, None
    if "qes:espresso" not in head and "<espresso" not in head:
        return False, None, None
    sch = re.search(r"(qes[-_][0-9A-Za-z_-]+)\.xsd", head)
    ver = re.search(r'<creator NAME="([^"]+)" VERSION="([^"]+)"', head)
    return True, (sch.group(1) if sch else None), (ver.group(2) if ver else None)


def _version_tuple(v):
    nums = re.findall(r"\d+", v or "")
    return tuple(int(n) for n in nums[:3]) if nums else ()


def _as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


_ELEMENTS = set("""H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn
Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu
Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm
Bk Cf Es Fm Md No Lr""".split())


def element_of(species_name):
    """Element symbol of a QE species label ('Fe', 'Fe1', 'Fe_up', 'fe2', 'O-' ...)."""
    m = re.match(r"[A-Za-z]+", species_name.strip())
    letters = m.group(0) if m else species_name.strip()
    for n in (2, 1):
        cand = letters[:n].capitalize()
        if cand in _ELEMENTS:
            return cand
    return letters


def _hubbard(dftu):
    """{'projection': ..., 'U_eV': {'Co 3d': 6.2, ...}, 'J_eV'/'V_eV' when present} or None."""
    if not isinstance(dftu, dict):
        return None
    out = {"projection": dftu.get("U_projection_type")}
    for tag, key in (("Hubbard_U", "U_eV"), ("Hubbard_J0", "J0_eV"), ("Hubbard_alpha", "alpha_eV")):
        vals = {}
        for item in _as_list(dftu.get(tag)):
            if isinstance(item, dict):
                vals[f"{item.get('@specie')} {item.get('@label', '')}".strip()] = round(float(item["$"]) * HARTREE_EV, 4)
        if vals:
            out[key] = vals
    if dftu.get("Hubbard_V"):
        out["V_eV"] = "intersite V present (not tabulated)"
    return out


class PwXml:
    """Decoded pw.x XML with accessors for the quantities MIDb needs."""

    def __init__(self, path):
        self.path = path
        ok, declared, version = sniff(path)
        if not ok:
            raise QEXMLError(f"{path} is not a pw.x XML file")
        self.declared_schema = declared
        self.qe_version = version
        vt = _version_tuple(version)
        if vt and vt < MIN_VERSION:
            raise QEXMLError(f"{path}: QE {version} is older than the supported {'.'.join(map(str, MIN_VERSION))}")
        self.schema_path = schema_for(declared)
        schema = _schema(self.schema_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                decoded = schema.to_dict(path, validation="lax")
            except Exception as error:   # truncated/corrupted file
                raise QEXMLError(f"{path}: cannot decode XML ({type(error).__name__}: {error})")
        data, errors = decoded if isinstance(decoded, tuple) else (decoded, [])
        self.n_validation_errors = len(errors)
        self.root = data
        if "output" not in self.root:
            raise QEXMLError(f"{path}: no <output> section (interrupted or initialization-only run)")
        self.input = self.root.get("input")          # absent in files rewritten by hp.x
        self.output = self.root["output"]
        bs = self.output.get("band_structure") or {}
        if "ks_energies" not in bs or "total_energy" not in self.output:
            raise QEXMLError(f"{path}: no eigenvalues/total energy (interrupted or initialization-only run)")

    # ------------------------------------------------------------------ identity
    @property
    def schema_name(self):
        return os.path.basename(self.schema_path)

    @property
    def calculation(self):
        """scf, nscf, bands, relax, vc-relax, md, vc-md; None when the input section is absent."""
        if not self.input:
            return None
        return self.input.get("control_variables", {}).get("calculation")

    @property
    def prefix(self):
        return (self.input or {}).get("control_variables", {}).get("prefix")

    @property
    def has_bands(self):
        return "band_structure" in self.output and "ks_energies" in self.output["band_structure"]

    # ------------------------------------------------------------------ structure
    def species_map(self):
        """{species label: element}, from the output (or input) atomic_species."""
        sp = (self.output.get("atomic_species") or (self.input or {}).get("atomic_species") or {})
        return {s["@name"]: element_of(s["@name"]) for s in _as_list(sp.get("species"))}

    def structure(self):
        """Final structure: labels, elements, cell (Å), cartesian positions (Å), alat (bohr)."""
        st = self.output["atomic_structure"]
        cell = np.array([st["cell"][f"a{i}"] for i in (1, 2, 3)], dtype=float) * BOHR_ANG
        atoms = _as_list(st["atomic_positions"]["atom"])
        labels = [a["@name"] for a in atoms]
        smap = self.species_map()
        return {"labels": labels, "elements": [smap.get(l, element_of(l)) for l in labels],
                "cell": cell, "positions": np.array([a["$"] for a in atoms], dtype=float) * BOHR_ANG,
                "alat_bohr": float(st.get("@alat", 0.0) or 0.0)}

    # ------------------------------------------------------------------ method
    def method(self):
        """Numerical parameters from the input section (Ry for cutoffs and degauss)."""
        inp = self.input or {}
        out = self.output
        basis = inp.get("basis") or {}
        obasis = out.get("basis_set") or {}
        ecutwfc = basis.get("ecutwfc", obasis.get("ecutwfc"))
        ecutrho = basis.get("ecutrho", obasis.get("ecutrho"))
        spin = inp.get("spin") or {}
        bs = out.get("band_structure", {})
        lsda = bool(spin.get("lsda", bs.get("lsda", False)))
        noncolin = bool(spin.get("noncolin", bs.get("noncolin", False)))
        spinorbit = bool(spin.get("spinorbit", bs.get("spinorbit", False)))
        bands = inp.get("bands") or {}
        occ = bands.get("occupations")
        if isinstance(occ, dict):
            occ = occ.get("$")
        sm = bands.get("smearing")
        smearing = None
        if occ == "smearing" and sm is not None:
            smearing = {"occupations": "smearing",
                        "smearing": sm.get("$") if isinstance(sm, dict) else sm,
                        "degauss": round(2 * float(sm.get("@degauss", 0.0)), 8) if isinstance(sm, dict) else None}
        kg = None
        kp = inp.get("k_points_IBZ") or {}
        mp = kp.get("monkhorst_pack")
        if isinstance(mp, dict):
            kg = [int(mp.get(f"@nk{i}", 1)) for i in (1, 2, 3)]
        species = _as_list((inp.get("atomic_species") or out.get("atomic_species") or {}).get("species"))
        dft = inp.get("dft") or out.get("dft") or {}
        xc = dft.get("functional")
        hubbard = _hubbard(dft.get("dftU") if isinstance(dft, dict) else None)
        return {
            "ecutwfc_Ry": round(2 * float(ecutwfc), 6) if ecutwfc is not None else None,
            "ecutrho_Ry": round(2 * float(ecutrho), 6) if ecutrho is not None else None,
            "xc": xc, "hubbard": hubbard,
            "pseudopotentials": {s["@name"]: s.get("pseudo_file") for s in species},
            "nspin": 2 if lsda else 1, "noncolin": noncolin, "spinorbit": spinorbit,
            "occupations": occ, "smearing": smearing, "kgrid": kg,
            "gamma_only": bool(basis.get("gamma_only", False)),
            "tot_charge": float(bands.get("tot_charge", 0.0) or 0.0),
        }

    # ------------------------------------------------------------------ energies
    def summary(self):
        """Total energy, Fermi level / HOMO / LUMO, magnetization (eV, μB), nelec."""
        out = self.output
        bs = out["band_structure"]
        res = {"etot_eV": float(out["total_energy"]["etot"]) * HARTREE_EV}
        for tag, key in (("fermi_energy", "efermi_eV"), ("highestOccupiedLevel", "homo_eV"),
                         ("lowestUnoccupiedLevel", "lumo_eV")):
            if bs.get(tag) is not None:
                res[key] = float(bs[tag]) * HARTREE_EV
        two = bs.get("two_fermi_energies")
        if two is not None:
            up, dw = (float(x) * HARTREE_EV for x in two)
            res.update(efermi_up_eV=up, efermi_dw_eV=dw, efermi_eV=max(up, dw))
        mag = out.get("magnetization") or {}
        if mag.get("total") is not None:
            res["magnetization_muB"] = float(mag["total"])
        if mag.get("absolute") is not None:
            res["abs_magnetization_muB"] = float(mag["absolute"])
        res["nelec"] = float(bs["nelec"])
        return res

    def eigenvalues(self):
        """Eigenvalues (eV) [spin][k][band] and k-points in cartesian 2π/alat units.

        Spin channels: 2 for LSDA, 1 otherwise (non-collinear bands are spinors)."""
        bs = self.output["band_structure"]
        ks = _as_list(bs["ks_energies"])
        kpts = np.array([k["k_point"]["$"] for k in ks], dtype=float)
        weights = np.array([float(k["k_point"].get("@weight", 0.0)) for k in ks])
        eig = np.array([k["eigenvalues"]["$"] for k in ks], dtype=float) * HARTREE_EV
        if bs.get("lsda"):
            nb = int(bs.get("nbnd_up") or eig.shape[1] // 2)
            eig = np.stack([eig[:, :nb], eig[:, nb:]])
        else:
            eig = eig[None]
        return eig, kpts, weights

    def kpoints_crystal(self):
        """k-points of the band structure in crystal (fractional) coordinates."""
        _, kcart, _ = self.eigenvalues()
        st = self.structure()
        cell_alat = st["cell"] / (st["alat_bohr"] * BOHR_ANG)   # lattice vectors in alat units
        return kcart @ cell_alat.T                                 # k_frac_i = k_cart · a_i (2π/alat units)


def read(path):
    return PwXml(path)
