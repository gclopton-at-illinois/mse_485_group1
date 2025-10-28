# scripts/make_vacancies.py
#!/usr/bin/env python3
"""
Create CeO2-δ by removing random O atoms from a stoichiometric supercell.

Inputs
------
- Either '--from <path/to/data.ceo2_*nm.lmp>' OR '--L 12 12 20' to locate the base file.
- One of '--delta <δ>' (0..0.5) or '--fracO <f_O_removed>' (0..1); related by f_O_removed = δ/2.
- '--seed <int>' for reproducibility.
- Optional '--neutralize {none,global}' (default: global) to keep total charge ≈ 0
  by spreading the missing O charge uniformly over all Ce.

Outputs
-------
- LAMMPS data: inputs/structure/data.ceo2_<Lx>x<Ly>x<Lz>nm_delta_<δ>.lmp
- JSON manifest with a 'vacancies' section recording δ, seed, and the removed O indices (1-based).
"""

from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
STRUCT_DIR = ROOT / "inputs" / "structure"

def fmt_dim(x: float) -> str:
    s = f"{x:.3f}".rstrip("0").rstrip(".")
    return (str(int(round(x))) if abs(x-round(x))<1e-9 else s).replace(".", "p")

def fmt_val(x: float) -> str:
    s = f"{x:.4f}".rstrip("0").rstrip(".")
    return s.replace(".", "p") if "." in s else s

def guess_paths(L: list[float]) -> tuple[Path, Path]:
    lx, ly, lz = (fmt_dim(L[0]), fmt_dim(L[1]), fmt_dim(L[2]))
    stem = f"data.ceo2_{lx}x{ly}x{lz}nm"
    return STRUCT_DIR / f"{stem}.lmp", STRUCT_DIR / f"{stem}.json"

def parse_args():
    p = argparse.ArgumentParser()
    gsrc = p.add_mutually_exclusive_group(required=False)
    gsrc.add_argument("--from", dest="from_path", type=Path,
                      help="Path to stoichiometric LAMMPS data file.")
    gsrc.add_argument("--L", nargs=3, type=float, metavar=("LX","LY","LZ"),
                      help="Box lengths in nm to locate data.ceo2_*nm.[lmp|json].")
    gdelta = p.add_mutually_exclusive_group(required=True)
    gdelta.add_argument("--delta", type=float, help="O deficiency δ in CeO2-δ (0..0.5).")
    gdelta.add_argument("--fracO", type=float, help="Fraction of O removed (0..1).")
    p.add_argument("--seed", type=int, required=True, help="RNG seed.")
    p.add_argument("--neutralize", choices=["none","global"], default="global",
                   help="Charge handling after removing O (default: global).")
    return p.parse_args()

def main():
    args = parse_args()

    # Resolve base file and manifest
    if args.from_path:
        base_lmp = args.from_path.resolve()
        base_json = base_lmp.with_suffix(".json")
    else:
        base_lmp, base_json = guess_paths(list(map(float, args.L)))
    if not base_lmp.exists() or not base_json.exists():
        sys.exit(f"Base files not found:\n  {base_lmp}\n  {base_json}")

    with open(base_json) as f:
        manifest = json.load(f)

    # δ ↔ fraction of O removed
    if args.delta is not None:
        if not (0.0 <= args.delta <= 0.5):
            sys.exit("delta must be in [0, 0.5].")
        fracO = args.delta / 2.0
        delta = args.delta
    else:
        if not (0.0 <= args.fracO <= 1.0):
            sys.exit("fracO must be in [0, 1].")
        fracO = args.fracO
        delta = 2.0 * fracO

    # Load atoms and identify species
    atoms = read(base_lmp.as_posix(), format="lammps-data", atom_style="charge")
    syms = np.array(atoms.get_chemical_symbols())
    q    = np.array(atoms.get_initial_charges(), dtype=float)
    n_tot = len(atoms)
    o_idx = np.nonzero(syms == "O")[0]
    ce_idx= np.nonzero(syms == "Ce")[0]
    nO0, nCe = len(o_idx), len(ce_idx)

    # How many O to remove?
    n_remove = int(round(fracO * nO0))
    if n_remove == 0:
        removed_1based = []
    else:
        rng = np.random.default_rng(args.seed)
        choose = rng.choice(o_idx, size=n_remove, replace=False)
        removed_1based = (choose + 1).tolist()   # record before deletion
        # Delete in descending order to preserve indices
        for idx in sorted(choose.tolist(), reverse=True):
            del atoms[idx]

    # Optional global neutralization on Ce
    if args.neutralize == "global" and n_remove > 0:
        # q_O is the (negative) O partial charge from the original manifest
        q_O = manifest.get("charges_by_species", {}).get("O", -1.2)
        total_missing = q_O * n_remove  # negative value removed -> system becomes +|q_O|*n_remove
        dq_per_Ce = - total_missing / nCe
        # Update charges on current atoms
        syms2 = np.array(atoms.get_chemical_symbols())
        q2 = np.array(atoms.get_initial_charges(), dtype=float)
        q2[syms2 == "Ce"] += dq_per_Ce
        atoms.set_initial_charges(q2)

    # Build output paths
    Lnm = manifest["box_lengths_nm"]["after_replication"]
    lx, ly, lz = (fmt_dim(Lnm[0]), fmt_dim(Lnm[1]), fmt_dim(Lnm[2]))
    out_stem = f"data.ceo2_{lx}x{ly}x{lz}nm_delta_{fmt_val(delta)}"
    out_lmp  = STRUCT_DIR / f"{out_stem}.lmp"
    out_json = STRUCT_DIR / f"{out_stem}.json"

    # Write data
    write(out_lmp.as_posix(), atoms, format="lammps-data", atom_style="charge", units="metal")

    # Update manifest
    new_manifest = dict(manifest)  # shallow copy ok
    syms_new = np.array(atoms.get_chemical_symbols())
    new_manifest["n_atoms"] = int(len(atoms))
    new_manifest["composition_counts"]["O"] = int(np.count_nonzero(syms_new == "O"))
    new_manifest["composition_counts"]["Ce"] = int(np.count_nonzero(syms_new == "Ce"))
    new_manifest["output_data"] = str(out_lmp.relative_to(ROOT))
    new_manifest["vacancies"] = {
        "delta": float(delta),
        "frac_O_removed": float(fracO),
        "n_removed": int(n_remove),
        "seed": int(args.seed),
        "removed_O_indices_1based": removed_1based,
        "neutralization": args.neutralize,
    }

    with open(out_json, "w") as f:
        json.dump(new_manifest, f, indent=2)

    print(f"[ok] Wrote {out_lmp}  (delta={delta:.4f}, removed {n_remove}/{nO0} O)")
    print(f"[ok] Manifest: {out_json}")

if __name__ == "__main__":
    main()
