#!/usr/bin/env python3
"""
# Target a given δ (per CeO2 formula unit)
python scripts/generate_vacancies.py \
    --input inputs/structure/data.ceo2_12x12x24nm.lmp \
    --delta 0.10 \
    --seed 12345

"""

from __future__ import annotations
import argparse, json, sys
from pathlib import Path

import numpy as np
from ase.io import read, write

# Repo paths (same pattern as make_supercell.py)
ROOT = Path(__file__).resolve().parents[1]
STRUCT_DIR = ROOT / "inputs" / "structure"


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate oxygen vacancies in a CeO2 LAMMPS data file to build CeO_{2-δ}."
    )
    p.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to parent LAMMPS data file (stoichiometric CeO2 or existing CeO_{2-δ}).",
    )

    vac_group = p.add_mutually_exclusive_group(required=True)
    vac_group.add_argument(
        "--delta",
        type=float,
        help="Target oxygen deficiency δ in CeO_{2-δ}. Used as N_vac ≈ δ * N_Ce.",
    )
    vac_group.add_argument(
        "--n-vac",
        type=int,
        help="Explicit number of oxygen vacancies to create.",
    )
    vac_group.add_argument(
        "--vac-frac",
        type=float,
        help="Fraction of oxygen atoms to remove (N_vac ≈ vac_frac * N_O).",
    )

    p.add_argument(
        "--ce-type",
        type=int,
        default=1,
        help="LAMMPS atom-type index for Ce (default: 1).",
    )
    p.add_argument(
        "--o-type",
        type=int,
        default=2,
        help="LAMMPS atom-type index for O (default: 2).",
    )

    p.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for vacancy placement (for reproducibility).",
    )

    p.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Optional explicit output path. If not given, the script appends "
            "'_delta<δlabel>' before .lmp on the input filename."
        ),
    )

    return p.parse_args()


def fmt_delta_label(delta: float) -> str:
    """Pretty δ label for filenames: 0.1 -> '0p1', 0.125 -> '0p125'."""
    x = abs(delta)
    s = f"{x:.4f}".rstrip("0").rstrip(".")
    return s.replace(".", "p") if s else "0"


def compute_vacancies(N_ce: int, N_o: int, args) -> tuple[int, float]:
    """
    Decide how many vacancies to create and report the effective δ.

    Returns
    -------
    N_vac : int
        Number of O atoms to remove.
    delta_eff : float
        Effective δ = 2 - N_O_final / N_Ce after rounding.
    """
    if N_ce <= 0 or N_o <= 0:
        sys.exit(f"[error] Non-positive species counts: N_Ce={N_ce}, N_O={N_o}.")

    if args.delta is not None:
        # δ per CeO_{2-δ}; ideal CeO2 has N_O = 2 * N_Ce
        N_vac_float = args.delta * N_ce
        N_vac = int(round(N_vac_float))
    elif args.n_vac is not None:
        N_vac = int(args.n_vac)
    else:
        # vac-frac case
        if not (0.0 < args.vac_frac < 1.0):
            sys.exit("[error] --vac-frac must be in (0,1).")
        N_vac_float = args.vac_frac * N_o
        N_vac = int(round(N_vac_float))

    if N_vac < 0:
        sys.exit(f"[error] Computed negative N_vac={N_vac}. Check inputs.")
    if N_vac == 0:
        sys.stderr.write("[warn] N_vac = 0; structure will be unchanged.\n")

    if N_vac > N_o:
        sys.exit(
            f"[error] Requested N_vac={N_vac} > N_O={N_o}. "
            "Cannot remove more oxygens than exist."
        )

    N_o_final = N_o - N_vac
    delta_eff = 2.0 - (N_o_final / float(N_ce))

    return N_vac, delta_eff


def choose_oxygen_sites_to_remove(o_indices: np.ndarray, N_vac: int, seed: int | None):
    rng = np.random.default_rng(seed)
    if N_vac == 0:
        return np.array([], dtype=int)
    return rng.choice(o_indices, size=N_vac, replace=False)


def main():
    args = parse_args()
    if not args.input.exists():
        sys.exit(f"[error] Input data file not found: {args.input}")

    # Read parent structure; ASE will guess H/He, but wee override that using types.
    atoms = read(args.input.as_posix(), format="lammps-data")

    # Extract LAMMPS 'type' array and force chemical symbols to Ce/O.
    try:
        types = atoms.get_array("type").astype(int)
    except Exception:
        sys.exit("[error] LAMMPS 'type' array not found in Atoms object.")

    symbols = []
    for t in types:
        if t == args.ce_type:
            symbols.append("Ce")
        elif t == args.o_type:
            symbols.append("O")
        else:
            sys.exit(
                f"[error] Encountered atom type {t} that is neither ce_type={args.ce_type} "
                f"nor o_type={args.o_type}. Refusing to guess."
            )
    atoms.set_chemical_symbols(symbols)

    # Now chemical symbols are Ce/O and we can count directly.
    syms = np.array(atoms.get_chemical_symbols())
    ce_mask = syms == "Ce"
    o_mask = syms == "O"
    N_ce = int(ce_mask.sum())
    N_o = int(o_mask.sum())
    N_tot = len(atoms)

    if N_ce == 0 or N_o == 0:
        sys.exit(
            f"[error] After type→symbol mapping, N_Ce={N_ce}, N_O={N_o}. "
            "Something is inconsistent with ce_type / o_type."
        )

    # Decide how many vacancies to create and the effective δ after rounding.
    N_vac, delta_eff = compute_vacancies(N_ce, N_o, args)

    # Choose which O indices to remove.
    o_indices = np.nonzero(o_mask)[0]
    vac_indices = choose_oxygen_sites_to_remove(o_indices, N_vac, args.seed)
    vac_indices_set = set(int(i) for i in vac_indices)

    # Build new Atoms object with those oxygens removed.
    keep_indices = [i for i in range(N_tot) if i not in vac_indices_set]
    atoms_vac = atoms[keep_indices]

    # counts
    syms_vac = np.array(atoms_vac.get_chemical_symbols())
    N_ce_new = int(np.count_nonzero(syms_vac == "Ce"))
    N_o_new = int(np.count_nonzero(syms_vac == "O"))

    if N_ce_new != N_ce:
        sys.exit(
            f"[bug] Ce count changed unexpectedly: before={N_ce}, after={N_ce_new}."
        )
    if N_o_new != N_o - N_vac:
        sys.exit(
            f"[bug] O count mismatch: before={N_o}, requested N_vac={N_vac}, "
            f"after={N_o_new}."
        )

    # Determine output path.
    if args.output is not None:
        out_path = args.output
    else:
        delta_label = fmt_delta_label(delta_eff)
        stem = args.input.stem
        out_name = f"{stem}_delta{delta_label}.lmp"
        out_path = args.input.with_name(out_name)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Write LAMMPS datar
    # specorder enforces element→type mapping on write; now syms = Ce/O so this is safe.
    write(
        out_path.as_posix(),
        atoms_vac,
        format="lammps-data",
        atom_style="charge",
        specorder=["Ce", "O"],
    )

    # Manifest JSON
    manifest = {
        "parent_data": (
            str(args.input.relative_to(ROOT))
            if str(args.input).startswith(str(ROOT))
            else str(args.input)
        ),
        "output_data": (
            str(out_path.relative_to(ROOT))
            if str(out_path).startswith(str(ROOT))
            else str(out_path)
        ),
        "rng_seed": args.seed,
        "parent_counts": {
            "N_atoms": int(N_tot),
            "N_Ce": int(N_ce),
            "N_O": int(N_o),
        },
        "type_mapping": {
            "Ce": int(args.ce_type),
            "O": int(args.o_type),
        },
        "vacancy_spec": {
            "requested_delta": args.delta,
            "requested_n_vac": args.n_vac,
            "requested_vac_frac": args.vac_frac,
        },
        "vacancy_results": {
            "N_vac": int(N_vac),
            "removed_o_indices": sorted(int(i) for i in vac_indices),
            "N_atoms_new": int(len(atoms_vac)),
            "N_Ce_new": int(N_ce_new),
            "N_O_new": int(N_o_new),
            "delta_eff": float(delta_eff),
        },
        "note": (
            "Ce/O identification is based on LAMMPS atom types; "
            "chemical symbols reset to Ce/O; atom_style=charge, units=metal."
        ),
    }

    manifest_path = out_path.with_suffix(".json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    # Summary
    print(
        f"[ok] Parent: {args.input}  "
        f"(N={N_tot}, N_Ce={N_ce}, N_O={N_o}, ce_type={args.ce_type}, o_type={args.o_type})"
    )
    print(
        f"[ok] Vacancies: N_vac={N_vac}, seed={args.seed}, "
        f"effective δ = {delta_eff:.6f}  → CeO_{{2-{delta_eff:.4f}}}"
    )
    print(
        f"[ok] New: {out_path}  "
        f"(N={len(atoms_vac)}, N_Ce={N_ce_new}, N_O={N_o_new})"
    )
    print(f"[ok] Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
