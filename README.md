# Quantum-Espresso-Parser-MIDb

`qe2midb` turns a folder with Quantum ESPRESSO outputs into `midb_qe.json`, the normalized
document read by the "Quantum ESPRESSO" entry type of [MIDb](https://midb.cloud)
(Materials Informatics Database). Structure, method, energetics, band structure, DOS,
projected DOS, phonon dispersion and DOS, Γ modes with IR intensities, ε∞ and Born charges,
for crystals and molecules. Format: [docs/midb_qe_format.md](docs/midb_qe_format.md).

Files are recognized **by content**, not by name or folder layout, and nothing comes from
outside the folder. Supported: QE **6.3 to 7.6**.

## Install and use

```bash
pip install quantum-espresso-parser-midb
```

```python
from qe2midb.normalize import normalize
doc, discovery = normalize("CALC_FOLDER")
```

Requirements: Python ≥ 3.9, numpy, xmlschema. Optional, installed separately: ASE (labels
of high-symmetry points; pymatgen is tried if ASE is missing) and spglib (space group). The official XML
schemas of pw.x (19 versions, from [QEF/qeschemas](https://github.com/QEF/qeschemas)) are
bundled in `qe2midb/schemas/`.

Example: `qe2midb examples/Si` (QE 7.6: scf, bands, nscf + dos.x + projwfc.x, ph.x +
q2r.x + matdyn.x + dynmat.x) gives `examples/Si/expected_midb_qe.json`.

## What is recognized

| file (any name, any sub-folder) | recognized by | used for |
|---|---|---|
| pw.x XML (`<prefix>.xml`, `<prefix>.save/data-file-schema.xml`) | root `qes:espresso`; type from `<calculation>` | scf: structure, method, energies; bands: band structure; nscf: band edges, DOS k-grid; relax/vc-relax: structure when there is no scf |
| XML rewritten by hp.x | no `<input>` section | ignored (reported) |
| dos.x output | header `E (eV) … EFermi` | total DOS |
| projwfc.x PDOS | `*.pdos_tot` + `*.pdos_atm#…` | PDOS per atom and per species/orbital |
| stdout of any QE program | `Program XXX v.Y starts` | versions; Löwdin spilling; IR intensities (dynmat.x) |
| inputs (pw, ph, matdyn, dynmat, dos, projwfc) | namelists | q-grid, `K_POINTS crystal_b` labels, DOS method, scf parameters |
| ph.x dynamical matrices | `Dynamical matrix file` | ε∞ and Born charges (E-U block) |
| matdyn `flfrq`, `fldos`; dynmat `filout` | `&plot nbnd` with 3·nat branches; `Frequency[cm^-1]`; single Γ block | phonon dispersion, phonon DOS, Γ modes |

- **XML**: read with the schema the file declares (nearest date when unpublished; QE 7.6
  declares `qes_260226`, available as `qes_dev_260226`), in lenient mode, because QE
  writes files with small deviations from its own schema.
- **k-path** (bands and phonons): linear axis in 1/Å rebuilt from the k-points; a
  discontinuity is a step more than 4× longer than both neighbouring steps (works for
  uniform-density paths and for `crystal_b` with a fixed number of points per segment);
  vertices at direction changes; labels from `crystal_b` comments, else ASE/pymatgen
  special points, else blank (the plot stays correct).
- **Crystal or molecule**: molecule when the k-grid is Γ only and there is more than 5 Å of
  vacuum in every direction; `--kind` overrides. The choice is recorded in `provenance`.
- **nscf/bands XML store a total energy of 0**: when the folder has only those, the energy
  comes from the scf stdout (a single converged run) or is left out, never 0; k-grid and
  smearing come from the scf input (`method.parameters_source`).

## Compatibility

XML values compared with the text output of the same run (total energy, Fermi level/HOMO,
number of atoms): |ΔE| ≤ 0.001 meV, |ΔE_F| < 0.05 meV.

| QE | declared schema | tested on | result |
|---|---|---|---|
| 6.1 | qes-20180510 | qeschema resources | not supported (older XML format) |
| 6.3 – 6.6 | 190304 – 200420 | aiida-quantumespresso fixtures | ok |
| 6.7MaX | 030920 | 2 production runs | ok |
| 6.8, 7.0 | 210716 | 52 production runs + fixtures | ok |
| 7.1 – 7.5 | 220603 – 250521 | 50 production runs + fixtures | ok |
| 7.6 | 260226 | MIDb demonstration set (16 systems) + fixtures | ok |

- Complete folders: 38 of 38 production folders (ZnO with and without DFT+U and Co/Fe
  doping, WO3 and CeO2 slabs, Pt nanoparticles, graphene…) and 50 of 86 aiida fixtures (the
  other 36 have no usable XML: failed runs with stdout only, truncated or
  initialization-only XML).
- PDOS: non-polarized, LSDA, non-collinear, spin-orbit (`p_j1.5`) and numbered species
  (`Fe1`); sum of projections equals `pdos_tot` within 0.2 %. k-resolved PDOS is skipped with
  a warning.
- On the 16 systems of the MIDb demonstration set the output equals that of the original
  campaign exporter.

## Tests

```bash
pip install ".[test]" && python -m pytest
```

Fixtures: pw.x XML of every schema from QE 6.3 to 7.6 and projwfc.x outputs from the
aiida-quantumespresso test suite (MIT, `tests/fixtures/qe2midb/aiida/LICENSE.txt`); matdyn
outputs from the Quantum ESPRESSO examples (GPL, see `tests/fixtures/qe2midb/qe_examples`);
`examples/Si` from the MIDb demonstration set. The Quantum ESPRESSO repository itself is not a
multi-version corpus: its example references and test-suite benchmarks are not regenerated
for each release.

## Known limitations

- QE < 6.3 (different XML format).
- One folder = one system; with several runs of the same kind, the most recent is used.
- Phonons only from the text outputs of matdyn.x / dynmat.x / ph.x (the ph.x XML is not used).
- Molecules without empty bands: HOMO-LUMO gap unknown; insulators without empty bands:
  reference at the VBM, gap unknown (warnings).
- Per-atom PDOS only up to 20 atoms by default (`--pdos-per-atom`); species/orbital sums are
  always kept.
- High-symmetry labels via pymatgen only when the cell is already the standard primitive one
  (ASE handles any cell).

## License

Copyright (C) 2026 James Almeida.

This program is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version. See [LICENSE](LICENSE).

Third-party files keep their own licenses: the aiida-quantumespresso fixtures in
`tests/fixtures/qe2midb/aiida/` (MIT, see the `LICENSE.txt` there) and the Quantum ESPRESSO
example outputs in `tests/fixtures/qe2midb/qe_examples/` (GPL-2.0-or-later). The XML schemas
in `qe2midb/schemas/` are distributed by the Quantum ESPRESSO Foundation
([QEF/qeschemas](https://github.com/QEF/qeschemas)), which does not state a license for them.

