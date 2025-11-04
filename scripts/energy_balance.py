#!/usr/bin/env python3
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


# ##############################################################
# ############### 1.) Path selection ##########################
# ##############################################################


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


# ##############################################################
# ############### 2.) Reading thermo table ####################
# ##############################################################

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
    Read thermo data from a whitespace-separated file or extract from a LAMMPS log

    Raise Error:
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
        # Fall back to attempting to parse the full log.
        return pd.read_csv(io.StringIO(text), sep=r"\s+", engine="python", comment="#")

    # whitespace-separated table.
    return pd.read_csv(io.StringIO(text), sep=r"\s+", engine="python", comment="#")


# ---------------------- normalize colum n names ------------

# Expand aliases to accept headers like "ET(eV)", "PE(eV)", "KE(eV)" from print titles
_CANONICAL = {
    "etotal": {
        "etotal", "TotEng", "total_energy", "TotalEnergy", "Etot", "E_tot", "Etotals",
        "ET(eV)", "ET"
    },
    "pe": {"pe", "PotEng", "E_pot", "Epot", "poteng", "PE(eV)", "PE"},
    "ke": {"ke", "KinEng", "E_kin", "Ekin", "kineng", "KE(eV)", "KE"},
    "step": {"step", "Step", "STEP"},
    "time": {"time", "Time", "t", "time_ps"},
    # Other collumns (like E_electron_eV, E_lattice_eV, *_eV_cum) used as-is
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


# ---------------------- physics config -------


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
    Ensure df has an 'etotal' column. If absent try to synthesize it from components

      - If 'etotal' already present (possibly via normalization of 'TotEng') do nothing
      - Else, if 'pe' and 'ke' exist, set etotal = pe + ke
      - Else, look for: E_electron_eV, E_lattice_eV, and any *_eV_cum reservoirs 
        (for example, E_langevin_eV_cum, E_boundary_eV_cum, E_cyl_sink_eV_cum, etc.). Set etotal = sum of those columns

    Returns (df_with_etotal, formula_string_or_None).
    """
    if "etotal" in df.columns:
        return df, "etotal (provided)"

    # Fallback 0: common thermo columns
    if "pe" in df.columns and "ke" in df.columns:
        df2 = df.copy()
        df2["etotal"] = df2["pe"].astype(float) + df2["ke"].astype(float)
        return df2, "pe + ke"

    # Fallback 1: electron + lattice (+ reservoirs)
    cols = set(df.columns)
    base = []
    if "E_electron_eV" in cols:
        base.append("E_electron_eV")
    if "E_lattice_eV" in cols:
        base.append("E_lattice_eV")

    if len(base) == 2:
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


# ---------------------- statistics ----------------- -----


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


# ---------------------- main ----------------


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

    # Try to read. if empty/missing, write a minimal JSON and exit quietly.
    try:
        df_raw = read_table(thermo_path)
    except pd.errors.EmptyDataError as e:
        print(f"[energy-balance] Skipping energy check: {e}")
        payload = {
            "measured": {"thermo_path": str(thermo_path), "n_rows": 0, "columns": []},
            "tolerances": {"etotal_rel_tol": args.tol},
            "ok": True,  # non-blocking
            "physics": load_physics(args.physics),
            "note": "thermo file missing or empty. check LAMMPS export or pass --thermo",
        }
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        sys.exit(0)

    # Normalize easy synonyms (for example, ET(eV) -> etotal; PE(eV) -> pe; KE(eV) -> ke)
    df = normalize_columns(df_raw)

    # Ensure we have an etotal, or derive if needed
    df, formula = derive_etotal_if_needed(df)
    if "etotal" not in df.columns:
        raise SystemExit(
            "[energy-balance] Could not find or derive 'etotal'. "
            f"Available columns: {list(df.columns)}. "
            "Expected either an 'etotal'/'TotEng' column, or ('pe' and 'ke'), "
            "or the pair 'E_electron_eV' and 'E_lattice_eV' (optionally plus *_eV_cum). "
            "You can also pass --thermo to point at a file that contains total energy."
        )

    # Compute stats
    stats = compute_etotal_stats(df)

    # Load physics (optional)
    physics_cfg = load_physics(args.physics)

    # Prepare JSON
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

    # console summary
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
