# midb_qe.json format

Document written by `qe2midb` and read by the "Quantum ESPRESSO" entry type of MIDb
(`MIDb/schema/qe/` in MIDb-FME).


Units: Å, eV, cm⁻¹ (THz where stated). **Energies of `bands`, `dos` and `pdos` are relative to
`energy_reference`**; absolute values are in `energetics`. Missing sections mean the step was not
run. Arrays are ordered `[spin][k][band]`; spin index 0 = up, 1 = down.

```
id (label, e.g. "Si"), formula (Hill), kind: "crystal" | "molecule"           [required]
method:                                                                          [required]
  code, qe_version, xml_schema, xc, pseudo_library, pseudopotentials {element: file},
  ecutwfc_Ry, ecutrho_Ry, nspin (1|2), kgrid_scf [3], kgrid_dos [3] | null, qgrid [3],
  smearing {occupations, smearing, degauss} | null, dos_method
structure:                                                                       [required]
  cell_angstrom [3][3], species [nat], positions_frac [nat][3],
  positions_cart_angstrom [nat][3], spacegroup (text), volume_angstrom3,
  lattice_abc_angstrom [3], angles_deg [3]
energetics (eV, absolute):                                                       [required]
  etot_eV, nelec, efermi_eV, homo_eV/lumo_eV (fixed occupations), vbm_eV/cbm_eV
  (gapped crystals), gap_eV (crystal: CBM-VBM or 0; molecule: LUMO-HOMO),
  magnetization_muB / abs_magnetization_muB (nspin = 2)
energy_reference: {name: "VBM" | "E_F" | "HOMO", value_eV}                       [required]
bands (crystals):
  path (e.g. "GXWKGLUWLK,UX"; comma = discontinuity), x [nk] (1/Å, linear axis;
  a repeated x marks a discontinuity), ticks [{x, label}] (two labels at the same x,
  e.g. K and U, at a discontinuity), kpoints_frac [nk][3], energies_eV [nspin][nk][nbnd]
dos:  energy_eV [ne], total [nspin][ne] (states/eV/cell), integrated [ne]
pdos:
  energy_eV [ne], entries [{atom_index, atom ("Ga1"), species, orbital (s|p|d|f),
  values [nspin][ne]}], species_orbital [{species, orbital, values [nspin][ne]}],
  lowdin_spilling
phonons:
  unit "cm-1", gamma_modes [{freq_cm1, freq_THz, ir_intensity ((D/Å)²/amu) | null,
  displacement [nat][3]}]
  crystals: path, x [nq], ticks, frequencies [nq][3·nat],
            dos {freq_cm1 [nf], total [nf], per_atom [nat][nf]}, min_freq_cm1,
            epsilon_inf [3][3] and born_charges [{atom, tensor [3][3]}] (insulators)
  molecules: n_vibrational_modes (3N−6 or 3N−5), vibrational_freqs_cm1 [n]
provenance: generated, run_dir, steps_done, qe_version_by_step, qe_version_expected
```

Notes on the derived document fields:
- `is_metal` is true for crystals referenced to E_F (e.g. PBE germanium, whose gap closes);
  `fermi_level` is set only in that case.
- Molecules get `lattice_type = "molecule"`, `dimensionality = 0` and their Schoenflies
  point group in `space_group`, so the landing-page aggregations stay consistent.
- Field names shared with other entry types (`formula`, `elements`, `space_group`,
  `final_energy`, `fermi_level`, `is_metal`, `is_spin_polarized`, ...) keep their meaning, so
  the search columns and the schema analysis work across types.

## Fields added by qe2midb (all optional for readers)

- `method`: `noncolin`, `spinorbit` (bool); `hubbard` (`null` or `{projection, U_eV: {"Co 3d": …}}`);
  `pseudo_library` may be `null` (not in the XML; set with `--pseudo-library`); `xml_schema` is
  the official schema used to read the file; `parameters_source` tells where the k-grid and
  smearing came from (`scf xml`, `nscf xml + scf input`, …).
- `energetics`: `efermi_up_eV` / `efermi_dw_eV` when the total magnetization is fixed;
  `gap_eV = null` without empty bands; `etot_eV` absent when only nscf/bands XML exist and no
  usable scf stdout (QE writes 0 in those XML files).
- `structure`: `spacegroup` from spglib (crystals; `null` without spglib and for molecules);
  `species_labels` when QE labels differ from elements (`Fe1`, `Fe2`).
- `pdos.entries` empty above 20 atoms by default; `species_orbital` always present.
- `provenance`: `generator`, `source_dir` (folder name only), `qe_version_by_program`,
  `xml {file, schema, declared_schema, validation_errors}`, `kind_source`
  (`heuristic`/`given`), `vacuum_angstrom`, `warnings`.
