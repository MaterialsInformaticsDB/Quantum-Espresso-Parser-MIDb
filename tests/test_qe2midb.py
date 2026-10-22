"""Tests of qe2midb, the generic QE -> midb_qe.json normalizer.

Cross-version fixtures: aiida-quantumespresso parser fixtures (MIT, see
fixtures/qe2midb/aiida/LICENSE.txt): pw.x XML of every schema from QE 6.3 to 7.6 with the
matching text output, and projwfc.x outputs (non-polarized, LSDA, non-collinear,
spin-orbit, numbered kinds, k-resolved). End-to-end regression: examples/Si (QE 7.6 run of
the MIDb demonstration set) against examples/Si/expected_midb_qe.json.

    python -m pytest
"""
import glob
import json
import os
import re
import shutil

import numpy as np
import pytest

from qe2midb import kpath as KP
from qe2midb import textread as T
from qe2midb.normalize import Discovery, hill_formula, normalize, vacuum_widths
from qe2midb.xmlread import HARTREE_EV, PwXml, QEXMLError, element_of, schema_for

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXAMPLE = os.path.join(ROOT, "examples", "Si")
EXPECTED = os.path.join(EXAMPLE, "expected_midb_qe.json")
AIIDA = os.path.join(HERE, "fixtures", "qe2midb", "aiida")
VERSION_DIRS = sorted(glob.glob(os.path.join(AIIDA, "pw", "default_xml_*")))
RY_EV = HARTREE_EV / 2


# ------------------------------------------------------------------ XML across versions
def test_schema_selection():
    assert os.path.basename(schema_for("qes_260226")) == "qes_dev_260226.xsd"   # QE 7.6
    assert os.path.basename(schema_for("qes_210716")) == "qes_210716.xsd"       # QE 6.8 / 7.0
    assert os.path.basename(schema_for("qes_241015")) == "qes_241104.xsd"       # nearest date


@pytest.mark.parametrize("d", VERSION_DIRS, ids=[os.path.basename(d) for d in VERSION_DIRS])
def test_every_schema_version_matches_text_output(d):
    """Value oracle: total energy, Fermi level/HOMO and nat from the XML vs pw.x stdout."""
    x = PwXml(os.path.join(d, "data-file-schema.xml"))
    s, st = x.summary(), x.structure()
    eig, k, w = x.eigenvalues()
    assert eig.ndim == 3 and eig.shape[1] == len(k) and np.isfinite(eig).all()
    out = os.path.join(d, "aiida.out")
    if not os.path.exists(out):
        pytest.skip("no stdout")
    t = open(out, errors="replace").read()
    e = re.findall(r"^!\s+total energy\s+=\s+(-?[\d.]+)\s+Ry", t, re.M)
    assert s["etot_eV"] == pytest.approx(float(e[-1]) * RY_EV, abs=1e-4)
    nat = int(re.search(r"number of atoms/cell\s*=\s*(\d+)", t).group(1))
    assert len(st["elements"]) == nat
    ef = re.findall(r"the Fermi energy is\s+(-?[\d.]+)\s+ev|highest occupied(?:, lowest unoccupied)? level(?:s)? \(ev\):\s+(-?[\d.]+)", t)
    if ef:
        level = s.get("efermi_eV", s.get("homo_eV"))
        assert level == pytest.approx(float(next(v for v in ef[-1] if v)), abs=1e-3)


def test_lsda_two_fermi_energies():
    x = PwXml(os.path.join(AIIDA, "pw", "tot_magnetization", "data-file-schema.xml"))
    s = x.summary()
    assert s["efermi_up_eV"] != s["efermi_dw_eV"] and s["efermi_eV"] == max(s["efermi_up_eV"], s["efermi_dw_eV"])
    assert x.eigenvalues()[0].shape[0] == 2 and x.method()["nspin"] == 2


@pytest.mark.parametrize("name", ["failed_interrupted_relax", "initialization_xml_new"])
def test_incomplete_xml_is_rejected(name):
    with pytest.raises(QEXMLError):
        PwXml(os.path.join(AIIDA, "pw", name, "data-file-schema.xml"))


def test_xml_without_input_is_ignored(tmp_path):
    """hp.x rewrites the XML without <input>: discovery must not take it for an scf."""
    src = open(os.path.join(VERSION_DIRS[-1], "data-file-schema.xml")).read()
    (tmp_path / "hp.xml").write_text(re.sub(r"<input>.*?</input>", "", src, flags=re.S))
    disc = Discovery(str(tmp_path))
    assert not disc.xml and "no <input>" in disc.ignored_xml[0][1]


def test_element_of_species_labels():
    assert [element_of(s) for s in ("Fe", "Fe1", "fe_up", "Co2", "O", "Cl1")] == ["Fe", "Fe", "Fe", "Co", "O", "Cl"]


# ------------------------------------------------------------------ projwfc variants
PDOS_CASES = {"nonpolarised": (1, False), "spinpolarised": (2, False), "noncollinear": (1, True),
              "spinorbit": (1, True), "numbered_kinds": (1, False)}


@pytest.mark.parametrize("case", sorted(PDOS_CASES))
def test_pdos_variants(case):
    channels, noncolin = PDOS_CASES[case]
    d = os.path.join(AIIDA, "projwfc", case)
    doc, disc = normalize(d)
    assert doc["method"]["noncolin"] is noncolin
    entries = doc["pdos"]["entries"]
    total = sum(np.array(e["values"]) for e in entries)
    assert total.shape[0] == channels
    assert {e["species"] for e in entries} <= {"Si", "Fe"}            # 'Si1' -> 'Si'
    assert {e["orbital"] for e in entries} <= {"s", "p", "d"}         # 'p_j1.5' -> 'p'
    ref = T.pdos_files(disc.pdos_prefixes[0], noncolin=noncolin)["pdos_sum"]
    if ref.max() > 0:
        assert np.abs(total - ref).max() / ref.max() < 0.01


def test_kresolved_pdos_is_skipped_with_warning():
    doc, _ = normalize(os.path.join(AIIDA, "projwfc", "kresolveddos"))
    assert "pdos" not in doc and any("k-resolved" in w for w in doc["provenance"]["warnings"])


# ------------------------------------------------------------------ helpers
def test_hill_formula():
    assert hill_formula(["H", "O", "H"]) == "H2O"
    assert hill_formula(["C"] * 6 + ["H"] * 6) == "C6H6"
    assert hill_formula(["N", "Ga", "N", "Ga"]) == "Ga2N2"


def test_vacuum_widths():
    cell = np.eye(3) * 16.0
    pos = [[8, 8, 8.3], [8, 8.77, 7.7], [8, 7.23, 7.7]]
    assert min(vacuum_widths(cell, pos)) > 14
    assert max(vacuum_widths(np.eye(3) * 3.0, [[0, 0, 0]])) == pytest.approx(3.0)


def test_kpath_reproduces_ase_path():
    """Axis, discontinuity (K|U) and labels rebuilt from the k-points alone."""
    d = json.load(open(EXPECTED))
    p = KP.build(d["bands"]["kpoints_frac"], d["structure"]["cell_angstrom"])
    assert np.allclose(p["x"], d["bands"]["x"], atol=1e-5)
    assert np.allclose([t["x"] for t in p["ticks"]], [t["x"] for t in d["bands"]["ticks"]], atol=1e-5)
    if not KP.special_points(d["structure"]["cell_angstrom"]):
        pytest.skip("labels need ASE or pymatgen")
    assert [t["label"] for t in p["ticks"]] == [t["label"] for t in d["bands"]["ticks"]]
    assert p["path"] == d["bands"]["path"]


def _crystal_b(sp, seq):
    """k-points as pw.x expands K_POINTS crystal_b: n points from each vertex to the next."""
    pts = []
    for (a, n), (b, _) in zip(seq[:-1], seq[1:]):
        pts += [sp[a] + (sp[b] - sp[a]) * i / n for i in range(n)]
    return np.array(pts + [sp[seq[-1][0]]])


@pytest.mark.parametrize("seq,path", [
    ([("G", 20), ("X", 20), ("W", 20), ("K", 20), ("G", 20), ("L", 1)], "GXWKGL"),
    ([("G", 20), ("X", 20), ("W", 20), ("K", 1), ("U", 20), ("X", 1)], "GXWK,UX"),      # weight 1 = jump
    ([("G", 40), ("X", 5), ("W", 40), ("L", 1)], "GXWL"),                               # uneven density
    ([("G", 3), ("X", 60), ("W", 1), ("L", 30), ("G", 1)], "GXW,LG"),
])
def test_kpath_crystal_b_paths(seq, path):
    cell = 5.43 / 2 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    sp = {"G": np.zeros(3), "X": np.array([0.5, 0, 0.5]), "W": np.array([0.5, 0.25, 0.75]),
          "K": np.array([0.375, 0.375, 0.75]), "L": np.array([0.5, 0.5, 0.5]), "U": np.array([0.625, 0.25, 0.625])}
    p = KP.build(_crystal_b(sp, seq), cell, labels=sp)                   # labels as from crystal_b comments
    assert p["path"] == path


def test_nscf_base_never_reports_zero_energy():
    """nscf XML stores etot = 0: without a usable scf stdout the energy must be absent."""
    doc, _ = normalize(os.path.join(AIIDA, "projwfc", "nonpolarised"))
    assert "etot_eV" not in doc["energetics"]
    assert any("total energy unknown" in w for w in doc["provenance"]["warnings"])


def test_matdyn_freq_reads_qe76_plot_format():
    q, w = T.matdyn_freq(os.path.join(EXAMPLE, "05_ph", "path.freq"))
    assert q.shape == (200, 3) and w.shape == (200, 6)
    assert w[0, 3:] == pytest.approx([504.2423] * 3)


# ------------------------------------------------------------------ regression on the campaign
ALLOWED = {"method.xml_schema", "method.noncolin", "method.spinorbit", "method.hubbard",
           "method.parameters_source", "provenance"}


def _diff(a, b, path=""):
    if path in ALLOWED:
        return []
    if isinstance(a, dict) and isinstance(b, dict):
        return [d for k in set(a) | set(b) for d in
                (_diff(a.get(k), b.get(k), f"{path}.{k}" if path else k) if k in a and k in b else [f"{path}.{k}"])
                if (f"{path}.{k}" if path else k) not in ALLOWED]
    if isinstance(a, list) and isinstance(b, list):
        try:
            A, B = np.array(a, dtype=float), np.array(b, dtype=float)
            return [] if A.shape == B.shape and (A.size == 0 or np.nanmax(np.abs(A - B)) < 2e-3) else [path]
        except (TypeError, ValueError):
            return [] if len(a) == len(b) and all(not _diff(x, y, path) for x, y in zip(a, b)) else [path]
    if isinstance(a, float) or isinstance(b, float):
        return [] if a is not None and b is not None and abs(a - b) < 2e-3 else [path]
    return [] if a == b else [path]


def test_example_end_to_end():
    """examples/Si (scf, bands, nscf/DOS/PDOS, DFPT phonons) -> the expected document."""
    doc, _ = normalize(EXAMPLE, system_id="Si", pseudo_library="SSSP PBE efficiency v2")
    assert _diff(json.load(open(EXPECTED)), doc) == []
    assert doc["energetics"]["gap_eV"] == pytest.approx(0.6107, abs=1e-4)
    assert doc["bands"]["path"] == "GXWKGLUWLK,UX" and doc["phonons"]["path"] == "GXWKGLUWLK,UX"
