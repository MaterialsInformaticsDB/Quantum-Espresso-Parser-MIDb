"""Build the normalized MIDb document (midb_qe.json) from a directory of QE outputs.

Files are recognized by content, not by name or folder layout: pw.x XML files are
classified by their <calculation>; dos.x, projwfc.x, ph.x, q2r.x, matdyn.x and dynmat.x
files by their headers or input namelists. Nothing comes from outside the directory.

    python -m qe2midb.normalize DIR [-o midb_qe.json] [--id NAME] [--kind crystal|molecule]
"""
import argparse
import datetime
import json
import os
import sys
from collections import defaultdict

import numpy as np

from qe2midb import kpath as KP
from qe2midb import textread as T
from qe2midb.xmlread import BOHR_ANG, PwXml, QEXMLError, sniff

FORMAT_VERSION = "1.0"
CM1_TO_THZ = 0.0299792458
MAX_TEXT_BYTES = 200 * 2**20
SKIP_DIR_SUFFIXES = (".save",)
VACUUM_MOLECULE = 5.0          # Å of vacuum in every direction -> molecule (with Gamma-only k)
PDOS_PER_ATOM_MAX = 20         # "auto": per-atom PDOS entries only up to this many atoms


def _r(a, nd=5):
    return np.round(np.asarray(a, dtype=float), nd).tolist()


# ====================================================================== discovery
class Discovery:
    """Inventory of the QE files of a directory, classified by content."""

    def __init__(self, directory, recursive=True):
        self.directory = os.path.abspath(directory)
        self.xml = defaultdict(list)          # calculation -> [PwXml]
        self.ignored_xml = []                 # (path, reason)
        self.text = defaultdict(list)         # kind -> [path]
        self.inputs = defaultdict(list)       # program -> [path]
        self.pdos_prefixes = []
        self.programs = {}                    # path -> (PROGRAM, version)
        self.warnings = []
        self._scan(recursive)

    def _files(self, recursive):
        for root, dirs, files in os.walk(self.directory):
            if root.endswith(SKIP_DIR_SUFFIXES):          # <prefix>.save: only its XML copy
                dirs[:] = []
                files = [f for f in files if f == "data-file-schema.xml"]
            elif not recursive:
                dirs[:] = [d for d in dirs if d.endswith(SKIP_DIR_SUFFIXES)]
            for f in files:
                yield os.path.join(root, f)

    def _scan(self, recursive):
        seen_xml = {}
        for p in sorted(self._files(recursive)):
            name = os.path.basename(p)
            try:
                size = os.path.getsize(p)
            except OSError:
                continue
            if name.endswith(".xml"):
                ok, _, _ = sniff(p)
                if not ok:
                    continue
                try:
                    doc = PwXml(p)
                except QEXMLError as error:
                    self.ignored_xml.append((p, str(error)))
                    continue
                if doc.calculation is None:
                    self.ignored_xml.append((p, "no <input> section (rewritten by hp.x or similar)"))
                    continue
                key = (doc.calculation, round(doc.summary()["etot_eV"], 6),
                       len(doc.output["band_structure"].get("ks_energies") or []))
                if key in seen_xml:                 # <prefix>.xml and data-file-schema.xml copies
                    continue
                seen_xml[key] = p
                self.xml[doc.calculation].append(doc)
                continue
            if size > MAX_TEXT_BYTES or "pdos_atm#" in name:
                continue
            if name.endswith(".pdos_tot"):
                self.pdos_prefixes.append(p[: -len(".pdos_tot")])
                continue
            h = T.head(p, 3000)
            if "\x00" in h:
                continue
            prog, ver = T.program_version(p)
            if prog:
                self.programs[p] = (prog, ver)
                self.text["stdout:" + prog].append(p)
                continue
            kind = T.qe_input_kind(p)
            if kind:
                self.inputs[kind].append(p)
            elif T.is_dos(p):
                self.text["dos"].append(p)
            elif T.is_dyn(p):
                self.text["dyn"].append(p)
            elif T.is_matdyn_freq(p):
                self.text["matdyn_freq"].append(p)
            elif T.is_phdos(p):
                self.text["phdos"].append(p)
            elif T.is_modes_file(p):
                self.text["modes"].append(p)
            elif name.endswith("dyn0") and T.dyn0_grid(p):
                self.text["dyn0"].append(p)

    # -------------------------------------------------------------- choices
    @staticmethod
    def latest(items, key=os.path.getmtime):
        return max(items, key=key) if items else None

    def pw(self, calculation):
        docs = self.xml.get(calculation, [])
        return max(docs, key=lambda d: os.path.getmtime(d.path)) if docs else None

    def report(self):
        return {"pw_xml": {k: [os.path.relpath(d.path, self.directory) for d in v] for k, v in self.xml.items()},
                "ignored_xml": [(os.path.relpath(p, self.directory), r) for p, r in self.ignored_xml],
                "text": {k: [os.path.relpath(p, self.directory) for p in v] for k, v in self.text.items()},
                "inputs": {k: [os.path.relpath(p, self.directory) for p in v] for k, v in self.inputs.items()},
                "pdos_prefixes": [os.path.relpath(p, self.directory) for p in self.pdos_prefixes],
                "warnings": self.warnings}


# ====================================================================== helpers
def hill_formula(elements):
    counts = defaultdict(int)
    for e in elements:
        counts[e] += 1
    order = (["C", "H"] + sorted(k for k in counts if k not in ("C", "H"))) if "C" in counts else sorted(counts)
    return "".join(f"{e}{counts[e] if counts[e] > 1 else ''}" for e in order if e in counts)


def vacuum_widths(cell, positions):
    """Empty space (Å) along each lattice direction, perpendicular to the other two vectors."""
    cell = np.asarray(cell, dtype=float)
    out = []
    for i in range(3):
        n = np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])
        n /= np.linalg.norm(n)
        width = abs(np.dot(cell[i], n))
        proj = np.sort((np.asarray(positions) @ n) % width)
        gaps = np.diff(np.concatenate([proj, [proj[0] + width]]))
        out.append(float(gaps.max()))
    return out


def guess_kind(method, cell, positions):
    k = method.get("kgrid")
    gamma = method.get("gamma_only") or k is None or list(k) == [1, 1, 1]
    vac = vacuum_widths(cell, positions)
    return ("molecule" if gamma and min(vac) > VACUUM_MOLECULE else "crystal"), vac


def band_edges(eig, nelec, noncolin=False):
    """VBM/CBM (eV) for integer filling without LSDA; None if metallic."""
    per_band = 1 if noncolin else 2
    nocc = int(round(nelec / per_band))
    if eig.shape[0] != 1 or abs(nelec / per_band - nocc) > 1e-6 or eig.shape[2] <= nocc:
        return None
    vbm, cbm = float(eig[0, :, nocc - 1].max()), float(eig[0, :, nocc].min())
    return (vbm, cbm) if cbm > vbm + 1e-3 else None


def is_linear(positions):
    if len(positions) <= 2:
        return True
    p = np.asarray(positions) - np.mean(positions, axis=0)
    return np.linalg.svd(p, compute_uv=False)[1] < 1e-3


def spacegroup(cell, frac, elements):
    try:
        import spglib
        nums = {e: i + 1 for i, e in enumerate(sorted(set(elements)))}
        ds = spglib.get_symmetry_dataset((cell, frac, [nums[e] for e in elements]), symprec=1e-3)
        if ds is None:
            return None
        symbol = ds.international if hasattr(ds, "international") else ds["international"]
        number = ds.number if hasattr(ds, "number") else ds["number"]
        return f"{symbol} ({number})"
    except Exception:  # noqa: BLE001 - optional
        return None


def frac_from_cart_alat(qcart, cell, alat_bohr):
    """q (cartesian, 2π/alat) -> crystal coordinates."""
    return np.asarray(qcart) @ (np.asarray(cell) / (alat_bohr * BOHR_ANG)).T


# ====================================================================== normalize
def normalize(directory, system_id=None, kind=None, pseudo_library=None, recursive=True, pdos_per_atom="auto"):
    """pdos_per_atom: "auto" (per-atom entries up to PDOS_PER_ATOM_MAX atoms), "all" or "none"."""
    disc = Discovery(directory, recursive)
    warn = disc.warnings
    base = disc.pw("scf") or disc.pw("vc-relax") or disc.pw("relax") or disc.pw("nscf") or disc.pw("bands")
    if base is None:
        why = "; ".join(f"{os.path.relpath(p, directory)}: {r}" for p, r in disc.ignored_xml[:3])
        raise ValueError(f"no usable pw.x XML in {directory}" + (f" ({why})" if why else ""))
    if disc.pw("scf") is None:
        warn.append(f"no scf XML; structure and energies from the {base.calculation} run")
    st, method, summ = base.structure(), base.method(), base.summary()
    method_source = f"{base.calculation} xml"
    if base.calculation in ("nscf", "bands"):
        # nscf/bands XML store etot = 0 and their own k-points/occupations: take the
        # energy from the scf stdout and the scf parameters from the scf input, if present
        summ.pop("etot_eV", None)
        energies = [e for e in (T.pw_stdout_scf_energy(p) for p, (prog, _) in disc.programs.items()
                                if prog == "PWSCF") if e is not None]
        if len(energies) == 1:
            summ["etot_eV"] = energies[0]
            warn.append("total energy taken from the scf stdout (nscf/bands XML store 0)")
        else:
            warn.append("total energy unknown: no single converged scf stdout in the directory")
        scf_inputs = [T.pw_input(p) for p in disc.inputs.get("pw", [])]
        scf_inputs = [i for i in scf_inputs if i["calculation"] == "scf"]
        if len(scf_inputs) == 1:
            method.update(kgrid=scf_inputs[0]["kgrid"], occupations=scf_inputs[0]["occupations"],
                          smearing=scf_inputs[0].get("smearing"))
            method_source = f"{base.calculation} xml + scf input"
        else:
            method.update(kgrid=None)
            warn.append("scf k-grid unknown (no single scf input in the directory)")
    cell, pos, elements = st["cell"], st["positions"], st["elements"]
    frac = pos @ np.linalg.inv(cell)
    kind_guess, vacuum = guess_kind(method, cell, pos)
    kind_source = "given" if kind else "heuristic"
    kind = kind or kind_guess
    if kind != kind_guess:
        warn.append(f"kind {kind!r} given, heuristic says {kind_guess!r} (vacuum {np.round(vacuum, 1).tolist()} Å)")
    crystal = kind == "crystal"
    smap = base.species_map()

    # ---- programs and versions
    versions = {}
    for p, (prog, ver) in disc.programs.items():
        versions.setdefault(prog, set()).add(ver)
    versions.setdefault("PWSCF(xml)", set()).update(d.qe_version for v in disc.xml.values() for d in v)
    versions = {k: sorted(v) for k, v in versions.items()}
    if len({v for vs in versions.values() for v in vs}) > 1:
        warn.append(f"more than one QE version in the directory: {versions}")

    nscf = disc.pw("nscf")
    # ---- phonon inputs/outputs
    ph_in = disc.latest(disc.inputs.get("ph", []))
    phin = T.ph_input(ph_in) if ph_in else None
    dyn0 = [p for p in disc.text.get("dyn0", [])]
    qgrid = (phin or {}).get("qgrid") or (T.dyn0_grid(dyn0[0]) if dyn0 else None)

    pseudos = {smap.get(k, k) if list(smap.values()).count(smap.get(k, k)) == 1 else k: v
               for k, v in method["pseudopotentials"].items()}
    doc = {
        "id": system_id or os.path.basename(os.path.abspath(directory)) or base.prefix,
        "formula": hill_formula(elements), "kind": kind,
        "method": {
            "code": f"Quantum ESPRESSO {base.qe_version}", "qe_version": base.qe_version,
            "xml_schema": base.schema_name, "xc": method["xc"],
            "pseudopotentials": pseudos, "pseudo_library": pseudo_library,
            "ecutwfc_Ry": method["ecutwfc_Ry"], "ecutrho_Ry": method["ecutrho_Ry"],
            "nspin": method["nspin"], "noncolin": method["noncolin"], "spinorbit": method["spinorbit"],
            "hubbard": method["hubbard"],
            "kgrid_scf": method["kgrid"],
            "kgrid_dos": (nscf.method()["kgrid"] if nscf is not None else None) if crystal else None,
            "qgrid": qgrid if crystal else ([1, 1, 1] if (phin or disc.text.get("modes")) else None),
            "smearing": method["smearing"],
            "dos_method": None, "parameters_source": method_source,
        },
        "structure": {
            "cell_angstrom": _r(cell, 6), "species": list(elements),
            "positions_frac": _r(frac, 6), "positions_cart_angstrom": _r(pos, 6),
            "spacegroup": spacegroup(cell, frac, elements) if crystal else None,
            "volume_angstrom3": round(abs(float(np.linalg.det(cell))), 4),
            "lattice_abc_angstrom": _r(np.linalg.norm(cell, axis=1), 5),
            "angles_deg": _r([np.degrees(np.arccos(np.dot(cell[j], cell[k]) / np.linalg.norm(cell[j]) / np.linalg.norm(cell[k])))
                              for j, k in ((1, 2), (0, 2), (0, 1))], 4),
        },
        "energetics": dict(summ),
    }
    if any(l != e for l, e in zip(st["labels"], elements)):
        doc["structure"]["species_labels"] = list(st["labels"])
    en = doc["energetics"]

    # ---- energy reference
    ref, ref_name = en.get("efermi_eV"), "E_F"
    if not crystal:
        ref, ref_name = en.get("homo_eV", en.get("efermi_eV")), "HOMO"
        if "lumo_eV" in en and "homo_eV" in en:
            en["gap_eV"] = en["lumo_eV"] - en["homo_eV"]
        else:
            en["gap_eV"] = None
            warn.append("molecule without empty bands (no LUMO): HOMO-LUMO gap unknown")

    # ---- bands
    bands = disc.pw("bands")
    eig_b = None
    if bands is not None and crystal:
        eig_b, _, _ = bands.eigenvalues()
    if crystal:
        sources = [d for d in (nscf, bands) if d is not None] or [base]
        edges = [band_edges(d.eigenvalues()[0], summ["nelec"], method["noncolin"]) for d in sources]
        if all(edges):
            en["vbm_eV"] = max(e[0] for e in edges)
            en["cbm_eV"] = min(e[1] for e in edges)
            en["gap_eV"] = en["cbm_eV"] - en["vbm_eV"]
            ref, ref_name = en["vbm_eV"], "VBM"
        elif method["occupations"] in ("fixed", None) and en.get("homo_eV") is not None:
            # insulator computed without empty bands: VBM known, gap unknown
            en["vbm_eV"] = en["homo_eV"]
            en["gap_eV"] = None
            ref, ref_name = en["vbm_eV"], "VBM"
            warn.append("fixed occupations without empty bands: band gap unknown")
        else:
            en["gap_eV"] = 0.0
    if ref is None:
        raise ValueError("no reference energy (Fermi level / HOMO) in the XML")
    if eig_b is not None:
        labels = {}
        for p in disc.inputs.get("pw", []):
            if "bands" in T.head(p, 20000).lower():
                labels.update({lab: k for k, lab in T.pw_kpoint_labels(p) if lab})
        path = KP.build(bands.kpoints_crystal(), cell, labels=labels or None)
        if path["n_labeled"] < path["n_vertices"]:
            warn.append(f"band path: {path['n_vertices'] - path['n_labeled']} of {path['n_vertices']} vertices without label")
        doc["bands"] = {"path": path["path"], "x": _r(path["x"], 6), "ticks": path["ticks"],
                        "kpoints_frac": _r(path["kpoints_frac"], 6), "energies_eV": _r(eig_b - ref, 4)}
    doc["energy_reference"] = {"name": ref_name, "value_eV": ref}

    # ---- DOS / PDOS
    dos_files = disc.text.get("dos", [])
    if dos_files:
        dd = T.dos_file(disc.latest(dos_files))
        doc["dos"] = {"energy_eV": _r(dd["energy"] - ref, 3), "total": _r(dd["total"], 5),
                      "integrated": _r(dd["integrated"], 4)}
    if disc.pdos_prefixes:
        prefix = disc.latest(disc.pdos_prefixes, key=lambda p: os.path.getmtime(p + ".pdos_tot"))
        try:
            pd = T.pdos_files(prefix, smap, noncolin=method["noncolin"])
            proj = disc.latest(disc.text.get("stdout:PROJWFC", []))
            keep_atoms = pdos_per_atom == "all" or (pdos_per_atom == "auto" and len(elements) <= PDOS_PER_ATOM_MAX)
            if not keep_atoms:
                warn.append(f"per-atom PDOS omitted ({len(elements)} atoms; species/orbital sums kept)")
            doc["pdos"] = {
                "energy_eV": _r(pd["energy"] - ref, 3),
                "entries": [{**{k: v for k, v in e.items() if k not in ("values", "label")},
                             "values": _r(e["values"], 5)} for e in pd["entries"]] if keep_atoms else [],
                "species_orbital": _sum_species(pd["entries"]),
                "lowdin_spilling": T.lowdin_spilling(proj) if proj else None,
            }
        except ValueError as error:
            warn.append(f"PDOS skipped: {error}")
    if nscf is not None and nscf.method()["occupations"] in ("tetrahedra", "tetrahedra_opt", "tetrahedra_lin"):
        doc["method"]["dos_method"] = nscf.method()["occupations"]
    else:
        for prog in ("dos", "projwfc"):
            path = disc.latest(disc.inputs.get(prog, []))
            nl = (T.namelist(open(path, errors="replace").read(), prog) or {}) if path else {}
            if nl.get("degauss"):
                doc["method"]["dos_method"] = f"gaussian {float(nl['degauss'].lower().replace('d', 'e')):g} Ry"
                break

    # ---- phonons
    ph = _phonons(disc, doc, cell, pos, elements, base, crystal, warn)
    if ph:
        doc["phonons"] = ph

    steps = sorted(set(disc.xml) | ({"ph"} if ph else set()) | ({"dos"} if "dos" in doc else set()))
    doc["provenance"] = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "generator": f"qe2midb {FORMAT_VERSION}", "source_dir": os.path.basename(os.path.abspath(directory)),
        "steps_done": steps, "qe_version_by_program": versions,
        "xml": {"file": os.path.relpath(base.path, directory), "schema": base.schema_name,
                "declared_schema": base.declared_schema, "validation_errors": base.n_validation_errors},
        "kind_source": kind_source, "vacuum_angstrom": _r(vacuum, 2),
        "warnings": warn,
        "energies_raw_note": "bands/dos/pdos referenciados a energy_reference; valores absolutos em energetics",
    }
    return doc, disc


def _sum_species(entries):
    acc = {}
    for e in entries:
        k = (e["species"], e["orbital"])
        acc[k] = acc[k] + e["values"] if k in acc else np.array(e["values"], dtype=float)
    return [{"species": s, "orbital": o, "values": _r(v, 5)} for (s, o), v in sorted(acc.items())]


def _phonons(disc, doc, cell, pos, elements, base, crystal, warn):
    ph = {"unit": "cm-1"}
    # Gamma modes: dynmat.x filout (single q block at Gamma)
    gamma = None
    for p in disc.text.get("modes", []):
        blocks = T.modes_file(p)
        if len(blocks) == 1 and np.allclose(blocks[0]["q"], 0.0, atol=1e-8):
            gamma = blocks[0]["modes"]
            break
    dyn_gamma = next((p for p in disc.text.get("dyn", []) if T.is_gamma_dyn(p)), None)
    eps, zs = T.epsilon_zstar(dyn_gamma) if dyn_gamma else (None, [])
    if gamma:
        ir = {}
        if crystal and zs:
            out = disc.latest(disc.text.get("stdout:DYNMAT", []))
            ir = {r["mode"]: r["ir"] for r in T.dynmat_ir(out)} if out else {}
        ph["gamma_modes"] = [{"freq_cm1": round(m["freq_cm1"], 3), "freq_THz": round(m["freq_THz"], 5),
                              "ir_intensity": ir.get(i + 1), "displacement": _r(m["eigvec"], 5)}
                             for i, m in enumerate(gamma)]
    if crystal:
        # bands.x writes the same '&plot nbnd=..' format: a phonon file has 3*nat branches
        freq_files = [p for p in disc.text.get("matdyn_freq", []) if T.plot_file_nbnd(p) == 3 * len(elements)]
        dos_runs = {os.path.basename(T.matdyn_input(p)["flfrq"]) for p in disc.inputs.get("matdyn", [])
                    if T.matdyn_input(p)["dos"]}
        path_files = [p for p in freq_files if os.path.basename(p) not in dos_runs] or freq_files
        if path_files:
            qcart, freqs = T.matdyn_freq(disc.latest(path_files))
            labels = {}
            for p in disc.inputs.get("matdyn", []):
                mi = T.matdyn_input(p)
                if mi["q_in_cryst_coord"]:
                    labels.update({lab: q for q, _, lab in mi["qpoints"] if lab})
            path = KP.build(frac_from_cart_alat(qcart, cell, base.structure()["alat_bohr"]), cell, labels or None)
            ph.update(path=path["path"], x=_r(path["x"], 6), ticks=path["ticks"], frequencies=_r(freqs, 3))
            ph["min_freq_cm1"] = round(float(freqs.min()), 2)
        dos_files = disc.text.get("phdos", [])
        if dos_files:
            d = T.phdos(disc.latest(dos_files), nat=len(elements))
            ph["dos"] = {"freq_cm1": _r(d["freq"], 2), "total": _r(d["total"], 6)}
            if d["per_atom"] is not None:
                ph["dos"]["per_atom"] = _r(d["per_atom"], 6)
        if eps:
            ph["epsilon_inf"] = _r(eps, 5)
            ph["born_charges"] = [{"atom": f"{s}{i + 1}", "tensor": _r(z, 5)}
                                  for i, (s, z) in enumerate(zip(elements, zs))]
    elif gamma:
        nzero = 5 if is_linear(pos) else 6
        freqs = sorted(m["freq_cm1"] for m in gamma)
        ph["n_vibrational_modes"] = len(elements) * 3 - nzero
        ph["vibrational_freqs_cm1"] = [round(f, 2) for f in freqs[nzero:]]
    return ph if len(ph) > 1 else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("directory")
    ap.add_argument("-o", "--output", help="output file (default: DIR/midb_qe.json)")
    ap.add_argument("--id", help="system label (default: name of the directory)")
    ap.add_argument("--kind", choices=("crystal", "molecule"), help="override the crystal/molecule heuristic")
    ap.add_argument("--pseudo-library", help="name of the pseudopotential library (metadata only)")
    ap.add_argument("--no-recursive", action="store_true", help="do not descend into sub-directories")
    ap.add_argument("--pdos-per-atom", choices=("auto", "all", "none"), default="auto",
                    help=f"per-atom PDOS entries (auto: only up to {PDOS_PER_ATOM_MAX} atoms)")
    ap.add_argument("--report", action="store_true", help="print the file inventory")
    a = ap.parse_args(argv)
    doc, disc = normalize(a.directory, a.id, a.kind, a.pseudo_library, recursive=not a.no_recursive,
                          pdos_per_atom=a.pdos_per_atom)
    out = a.output or os.path.join(a.directory, "midb_qe.json")
    with open(out, "w") as f:
        json.dump(doc, f, ensure_ascii=False, separators=(",", ":"))
    sections = [k for k in ("bands", "dos", "pdos", "phonons") if k in doc]
    print(f"{doc['id']}: {doc['formula']} {doc['kind']} QE {doc['method']['qe_version']} -> {out} "
          f"({os.path.getsize(out) / 1024:.0f} KiB; {', '.join(sections) or 'structure/energies only'})")
    for w in doc["provenance"]["warnings"]:
        print("  aviso:", w)
    if a.report:
        print(json.dumps(disc.report(), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
