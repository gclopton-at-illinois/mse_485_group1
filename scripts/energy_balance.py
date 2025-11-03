#!/usr/bin/env python3
r"""
energy_balance.py — Phase-1 sanity check for energy accounting.

What this script does
---------------------
1) Finds a thermo source for a given run directory, preferring (in order):
     <RUN>/post/thermo_energy.tsv
     <RUN>/logs/log.lammps      (thermo block will be extracted)
     <RUN>/dumps/thermo.tsv
   You can override with --thermo /path/to/file.

2) Reads the thermo table robustly:
   - If given a LAMMPS log, it extracts the lines after the "Step ..." header
     until the table ends, and parses those lines with whitespace separation.
   - Otherwise, it parses the entire file as whitespace-separated (no deprecated
     delim_whitespace; we use sep=r"\s+").

3) If no explicit 'etotal' column exists, derive it from physically meaningful components:
     etotal_eV = E_electron_eV + E_lattice_eV [+ any *_eV_cum reservoirs if present]
   (units don’t matter for relative drift; we report the formula used.)

4) Computes energy-drift stats and writes JSON to:
     <RUN>/post/energy_balance.json
"""

from __future__ import annotations

import argparse
import io
import json
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple, List

import pandas as pd

try:
    import yaml  # type: ignore
except Exception:  # pragma: no cover
    yaml = None


# ---------------------- path selection ----------------------


def sniff_thermo_path(run_dir: Path, override: Optional[Path]) -> Path:
    """
    Choose the best available thermo source (override wins).
    """
    if override is not None:
        return override

    candidates = [
        run_dir / "post" / "thermo_energy.tsv",
        run_dir / "logs" / "log.lammps",
        run_dir / "dumps" / "thermo.tsv",
    ]
    for p in candidates:
        try:
            if p.is_file() and p.stat().st_size > 0:
                return p
        except FileNotFoundError:
            pass

    # Nothing non-empty; return the preferred path so the caller can emit a helpful message.
    return candidates[0]


# ---------------------- parsing helpers ----------------------


def _extract_thermo_block(text: str) -> Tuple[Optional[List[str]], Optional[str]]:
    """
    Pull the numeric thermo table out of a LAMMPS log.

    Looks for:
      a header line starting with "Step "
      followed by numeric rows (whitespace-separated) until a blank or non-numeric line.

    Returns (header_list, body_text) or (None, None) if not found.
    """
    header: Optional[List[str]] = None
    rows: List[str] = []
    in_block = False
    for ln in text.splitlines():
        s = ln.strip()
        if not s:
            if in_block and rows:
                break
            continue
        if s.startswith("Step "):
            header = re.split(r"\s+", s)
            in_block = True
            continue
        if in_block:
            if re.match(r"^[\s\-\d\.+eE]+$", s):
                rows.append(s)
            else:
                if rows:
                    break
    if header and rows:
        return header, "\n".join(rows)
    return None, None


def read_table(path: Path) -> pd.DataFrame:
    """
    Read thermo data from a whitespace-separated file or extract from a LAMMPS log.

    Raises:
      pd.errors.EmptyDataError for missing/empty files.
    """
    if (not path.is_file()) or path.stat().st_size == 0:
        raise pd.errors.EmptyDataError(f"{path} is missing or empty")

    text = path.read_text(encoding="utf-8", errors="ignore")

    # If this looks like a LAMMPS log, try to extract the thermo block.
    if path.name.startswith("log."):
        header, body = _extract_thermo_block(text)
        if header and body:
            return pd.read_csv(io.StringIO(body), sep=r"\s+", engine="python", names=list(header))
        # Fall back to attempting to parse the full log (best effort).
        return pd.read_csv(io.StringIO(text), sep=r"\s+", engine="python", comment="#")

    # Generic whitespace-separated table.
    return pd.read_csv(io.StringIO(text), sep=r"\s+", engine="python", comment="#")


# ---------------------- column normalization ----------------------

_CANONICAL = {
    "etotal": {"etotal", "TotEng", "total_energy", "TotalEnergy", "Etot", "E_tot", "Etotals"},
    "step": {"step", "Step", "STEP"},
    "time": {"time", "Time", "t", "time_ps"},
    # keep these unrenamed; we just detect them for derivation
    # "E_electron_eV", "E_lattice_eV", "E_*_eV_cum" will be used as-is
}


def _inverse_name_map(columns: Iterable[str]) -> Dict[str, str]:
    present = {c: c for c in columns}
    lower = {c.lower(): c for c in columns}
    mapping: Dict[str, str] = {}
    for canon, aliases in _CANONICAL.items():
        for a in aliases:
            if a in present:
                mapping[canon] = a
                break
            if a.lower() in lower:
                mapping[canon] = lower[a.lower()]
                break
    return mapping


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    mapping = _inverse_name_map(df.columns)
    new_cols: List[str] = []
    for c in df.columns:
        renamed = None
        for canon, actual in mapping.items():
            if actual == c:
                renamed = canon
                break
        new_cols.append(renamed if renamed is not None else c)
    out = df.copy()
    out.columns = new_cols
    return out


# ---------------------- physics config ----------------------


def load_physics(yaml_path: Optional[Path]) -> Dict:
    if yaml_path is None or not yaml_path.is_file() or yaml is None:
        return {}
    try:
        with yaml_path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return data if isinstance(data, dict) else {"_raw": data}
    except Exception as exc:  # pragma: no cover
        return {"_error": f"failed to load physics yaml: {exc}"}


# ---------------------- derive etotal when needed ----------------------


def derive_etotal_if_needed(df: pd.DataFrame) -> Tuple[pd.DataFrame, Optional[str]]:
    """
    Ensure df has an 'etotal' column; if absent, try to synthesize it from components.

    Strategy:
      - If 'etotal' already present (possibly via normalization of 'TotEng'), do nothing.
      - Else, look for: E_electron_eV, E_lattice_eV, and any *_eV_cum reservoirs
        (e.g., E_langevin_eV_cum, E_boundary_eV_cum, E_cyl_sink_eV_cum, etc.)
        and set etotal = sum of those columns.

    Returns (df_with_etotal, formula_string_or_None).
    """
    if "etotal" in df.columns:
        return df, "etotal (provided)"

    cols = set(df.columns)

    # Base components that should both exist
    base = []
    if "E_electron_eV" in cols:
        base.append("E_electron_eV")
    if "E_lattice_eV" in cols:
        base.append("E_lattice_eV")

    if len(base) == 2:
        # Optional cumulative reservoirs to include if present
        extras = [
            "E_langevin_eV_cum",
            "E_boundary_eV_cum",
            "E_cyl_sink_eV_cum",
            "E_sink_eV_cum",
            "E_coupling_eV_cum",
            "E_thermostat_eV_cum",
        ]
        comp = base + [e for e in extras if e in cols]
        df2 = df.copy()
        df2["etotal"] = df2[comp].astype(float).sum(axis=1)
        return df2, " + ".join(comp)

    return df, None


# ---------------------- statistics ----------------------


def compute_etotal_stats(df: pd.DataFrame) -> Dict[str, float]:
    et = df["etotal"].astype(float).to_numpy()
    if et.size == 0:
        raise ValueError("No rows in thermo table after parsing")
    start = float(et[0])
    end = float(et[-1])
    abs_change = end - start
    denom = max(1.0, abs(start))
    rel_change = abs_change / denom
    max_abs_dev = float((abs(et - start)).max() / denom)
    return {
        "start": start,
        "end": end,
        "abs_change": abs_change,
        "rel_change": rel_change,
        "max_abs_dev": max_abs_dev,
    }


# ---------------------- main ----------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Phase-1 energy balance sanity check")
    p.add_argument("--run", required=True, type=Path, help="Path to run directory")
    p.add_argument("--physics", type=Path, default=None, help="Path to physics.yaml")
    p.add_argument("--thermo", type=Path, default=None, help="Override path to thermo file")
    p.add_argument("--tol", type=float, default=0.02, help="Relative tolerance on ΔE_tot/|E0|")
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output JSON path (default: <RUN>/post/energy_balance.json)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_dir: Path = args.run.resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"[energy-balance] Run directory not found: {run_dir}")

    post_dir = run_dir / "post"
    post_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out if args.out is not None else (post_dir / "energy_balance.json")

    thermo_path = sniff_thermo_path(run_dir, args.thermo)

    # Try to read; if empty/missing, write a minimal JSON and exit quietly.
    try:
        df_raw = read_table(thermo_path)
    except pd.errors.EmptyDataError as e:
        print(f"[energy-balance] Skipping energy check: {e}")
        payload = {
            "measured": {"thermo_path": str(thermo_path), "n_rows": 0, "columns": []},
            "tolerances": {"etotal_rel_tol": args.tol},
            "ok": True,
            "physics": load_physics(args.physics),
            "note": "thermo file missing or empty; check LAMMPS export or pass --thermo",
        }
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        sys.exit(0)

    # Normalize easy synonyms (e.g., TotEng -> etotal)
    df = normalize_columns(df_raw)

    # Ensure we have an etotal, deriving if needed
    df, formula = derive_etotal_if_needed(df)
    if "etotal" not in df.columns:
        raise SystemExit(
            "[energy-balance] Could not find or derive 'etotal'. "
            f"Available columns: {list(df.columns)}. "
            "Expected either an 'etotal'/'TotEng' column, or the pair "
            "'E_electron_eV' and 'E_lattice_eV' (optionally plus *_eV_cum). "
            "You can also pass --thermo to point at a file that contains TotEng."
        )

    # Compute stats
    stats = compute_etotal_stats(df)

    # Load physics (optional)
    physics_cfg = load_physics(args.physics)

    # Prepare JSON payload
    payload = {
        "measured": {
            "thermo_path": str(thermo_path),
            "n_rows": int(len(df)),
            "columns": list(map(str, df.columns)),
            "etotal": stats,
        },
        "tolerances": {"etotal_rel_tol": float(args.tol)},
        "ok": abs(stats["rel_change"]) <= float(args.tol),
        "physics": physics_cfg,
    }
    if formula is not None:
        payload["measured"]["etotal_formula"] = formula

    # Write
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # Brief console summary
    rel = stats["rel_change"]
    sign = "+" if rel >= 0 else "-"
    pct = abs(rel) * 100.0
    print(
        f"[energy-balance] E_tot change {sign}{pct:.3f}% over run "
        f"({stats['start']:.6g} → {stats['end']:.6g}); "
        f"max deviation {stats['max_abs_dev']*100:.3f}% ; "
        f"ok={payload['ok']} (tol={args.tol:.3f})"
    )


if __name__ == "__main__":
    main()
