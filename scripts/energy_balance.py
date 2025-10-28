# scripts/energy_balance.py
#!/usr/bin/env python3
"""
energy_balance.py — Phase-1 sanity gate for source normalization.

What it does
------------
Reads a run directory, pulls the intended lineal energy deposition S_e from
`config/physics.yaml` (or a snapshot), and compares it to the *measured*
deposition inferred from time series:
  measured_per_nm ≈ Δ(E_electron + E_lattice) + E_langevin_cum  (all per nm)

It writes `post/energy_balance.json` with the numbers and prints a PASS/FAIL
line (tolerance default 2%).

Inputs it looks for (in order)
------------------------------
1) `--thermo` file if you pass one (CSV/TSV with headers)
2) `<run>/dumps/thermo.tsv`
The file must contain either:
  (A) a column named `E_dep_eV_per_nm_cum` (then we just take its final value), or
  (B) columns `E_electron_eV`, `E_lattice_eV`, and optionally `E_langevin_eV_cum`
      (cumulative thermostat work; sign may be ±). We compute:
          ΔE = (Ee+El)_final - (Ee+El)_initial
          measured_total_eV = ΔE + (E_langevin_eV_cum_final or 0)
      then divide by Lz[nm] from the structure manifest.

If the columns aren’t present, the script explains what to enable in thermo
output and exits with a non-zero status (useful in CI).

Usage
-----
python scripts/energy_balance.py --run runs/phase1_baseline/<stamp>_... \
                                 --physics config/physics.yaml \
                                 --tol 0.02

You can override the column names if your thermo headers differ:
  --col-ee E_electron_eV --col-el E_lattice_eV --col-rim E_langevin_eV_cum
"""

from __future__ import annotations
import argparse, json, sys, re
from pathlib import Path

import yaml
import pandas as pd

EV_TO_J = 1.602176634e-19

def find_structure_manifest(run_dir: Path) -> Path | None:
    # Prefer the snapshot inside the run
    for base in ("inputs_snapshot/structure", "inputs/structure"):
        p = run_dir / base
        if p.is_dir():
            cands = sorted(p.glob("*.json"))
            if cands:
                return cands[0]
    return None

def read_box_length_nm(manifest_path: Path) -> float:
    data = json.loads(manifest_path.read_text())
    try:
        return float(data["box_lengths_nm"]["after_replication"][2])
    except Exception:
        raise SystemExit(f"[energy_balance] Could not parse Lz from {manifest_path}")

def load_physics_yaml(physics_path: Path) -> dict:
    with open(physics_path, "r") as f:
        return yaml.safe_load(f)

def sniff_thermo_path(run_dir: Path, override: Path | None) -> Path:
    if override:
        return override
    p = run_dir / "dumps" / "thermo.tsv"
    if p.exists():
        return p
    raise SystemExit(f"[energy_balance] No thermo file found. Looked for {p}. "
                     "Pass --thermo <path> or enable writing dumps/thermo.tsv with the required columns.")

def read_table(path: Path) -> pd.DataFrame:
    # Try TSV first; fall back to any whitespace-separated
    try:
        df = pd.read_csv(path, sep=None, engine="python", comment="#")
    except Exception:
        df = pd.read_csv(path, delim_whitespace=True, comment="#")
    if df.empty:
        raise SystemExit(f"[energy_balance] {path} is empty.")
    # Normalize header whitespace
    df.columns = [re.sub(r"\s+", "_", str(c).strip()) for c in df.columns]
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, required=True, help="Run directory (contains dumps/, logs/, post/).")
    ap.add_argument("--physics", type=Path, default=Path("config/physics.yaml"),
                    help="Physics card (for S_e). Defaults to config/physics.yaml relative to CWD.")
    ap.add_argument("--thermo", type=Path, help="Override path to thermo/energy time series (CSV/TSV).")
    ap.add_argument("--tol", type=float, default=0.02, help="Relative tolerance for PASS (e.g., 0.02 for 2%).")
    ap.add_argument("--col-ee", default="E_electron_eV", help="Column name for electron energy [eV].")
    ap.add_argument("--col-el", default="E_lattice_eV", help="Column name for lattice energy [eV].")
    ap.add_argument("--col-rim", default="E_langevin_eV_cum", help="Column for cumulative rim work [eV].")
    ap.add_argument("--col-dep", default="E_dep_eV_per_nm_cum", help="Column for cumulative deposited energy per nm [eV/nm].")
    args = ap.parse_args()

    run_dir = args.run.resolve()
    post_dir = run_dir / "post"
    post_dir.mkdir(parents=True, exist_ok=True)

    # Expected per-nm from physics.yaml: Se [eV/Å] × 10 Å/nm
    phys = load_physics_yaml(args.physics.resolve())
    try:
        Se_eV_per_A = float(phys["source"]["Se_eV_per_A"])
    except Exception:
        raise SystemExit(f"[energy_balance] Could not find source.Se_eV_per_A in {args.physics}")
    expected_eV_per_nm = Se_eV_per_A * 10.0
    expected_J_per_m = expected_eV_per_nm * EV_TO_J * 1e9

    # Geometry
    man_path = find_structure_manifest(run_dir)
    if not man_path:
        raise SystemExit("[energy_balance] Could not find a structure manifest (…/inputs_snapshot/structure/*.json).")
    Lz_nm = read_box_length_nm(man_path)

    # Load thermo/energy table
    thermo_path = sniff_thermo_path(run_dir, args.thermo)
    df = read_table(thermo_path)

    # Compute measured per-nm
    measured_eV_per_nm = None
    if args.col_dep in df.columns:
        measured_eV_per_nm = float(df[args.col_dep].iloc[-1])
        basis = f"direct from column {args.col_dep}"
    else:
        if args.col-ee not in df.columns or args.col-el not in df.columns:
            raise SystemExit(f"[energy_balance] Missing required columns in {thermo_path}."
                             f" Need either '{args.col_dep}' OR both '{args.col-ee}' and '{args.col-el}'.")
        Ee0, Ee1 = float(df[args.col-ee].iloc[0]), float(df[args.col-ee].iloc[-1])
        El0, El1 = float(df[args.col-el].iloc[0]), float(df[args.col-el].iloc[-1])
        dE = (Ee1 + El1) - (Ee0 + El0)  # eV (total, full box)
        rim = float(df[args.col-rim].iloc[-1]) if args.col-rim in df.columns else 0.0
        measured_total_eV = dE + rim
        measured_eV_per_nm = measured_total_eV / Lz_nm
        basis = (f"Δ(Ee+El) + {args.col-rim if args.col-rim in df.columns else '0'}; "
                 f"divided by Lz={Lz_nm:.3f} nm")

    measured_J_per_m = measured_eV_per_nm * EV_TO_J * 1e9
    rel_err = abs(measured_eV_per_nm - expected_eV_per_nm) / max(expected_eV_per_nm, 1e-30)
    pass_bool = rel_err <= args.tol

    # Write JSON artifact
    out = {
        "expected": {
            "Se_eV_per_A": Se_eV_per_A,
            "per_nm_eV": expected_eV_per_nm,
            "per_m_J": expected_J_per_m,
        },
        "measured": {
            "per_nm_eV": measured_eV_per_nm,
            "per_m_J": measured_J_per_m,
            "basis": basis,
            "thermo_path": str(thermo_path),
        },
        "geometry": {
            "Lz_nm": Lz_nm,
            "structure_manifest": str(man_path),
        },
        "tolerance_rel": args.tol,
        "relative_error": rel_err,
        "pass": bool(pass_bool),
    }
    out_path = post_dir / "energy_balance.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    status = "PASS" if pass_bool else "FAIL"
    print(f"[energy_balance] {status}  "
          f"measured={measured_eV_per_nm:.3f} eV/nm  "
          f"expected={expected_eV_per_nm:.3f} eV/nm  "
          f"(rel.err={100*rel_err:.2f}% ; basis: {basis})")
    sys.exit(0 if pass_bool else 2)

if __name__ == "__main__":
    main()
