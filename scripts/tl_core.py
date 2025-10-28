#!/usr/bin/env python3
"""
Compute lattice core temperature T_l(t) from LAMMPS dumps (metal units),
restricting atoms to a cylinder r <= r_core_nm about the box center.

Inputs (required)
-----------------
--pattern runs/.../dumps/atoms.spike.*.lammpstrj
  Each frame must have columns: id type x y z vx vy vz

Outputs
-------
CSV at --out with columns: t_ps, Tl_core_K
Also writes a tiny meta file alongside the CSV with basic parameters.

Optional (future)
-----------------
If a sibling file 'Te_core.csv' appears later (e.g., from TTM grid averaging),
make_phase1_figs.py will automatically overlay Te and Tl. This script does not
attempt to parse the TTM grid itself (Phase-1 keeps atoms-only here).
"""

import argparse, os, glob, re
import numpy as np, pandas as pd

# ---- constants for LAMMPS metal units ----
KB_eV_per_K = 8.617333262e-5      # eV/K
MVV2E       = 1.0364269e-4        # eV per (amu * (Å/ps)^2)
A_PER_NM    = 10.0

# Type → mass (amu) map (matches input deck: mass 1 140.116; mass 2 15.999)
MASS_AMU = {1: 140.116, 2: 15.999}

def _step_key(path: str) -> int:
    m = re.search(r"(\d+)(?:\D*)$", os.path.basename(path))
    return int(m.group(1)) if m else 0

def parse_sequence(pattern):
    """Yield frames as (box, ids, types, pos, vel) for all files matching pattern in numeric order."""
    paths = sorted(glob.glob(pattern), key=_step_key)
    if not paths:
        raise SystemExit(f"No spike dumps matched: {pattern}")
    for path in paths:
        with open(path, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if not line.startswith("ITEM: TIMESTEP"):
                    continue
                _timestep = int(f.readline().strip())
                assert f.readline().startswith("ITEM: NUMBER OF ATOMS")
                natoms = int(f.readline().strip())
                bounds_hdr = f.readline().strip()
                assert bounds_hdr.startswith("ITEM: BOX BOUNDS")
                xlo,xhi = map(float, f.readline().split()[:2])
                ylo,yhi = map(float, f.readline().split()[:2])
                zlo,zhi = map(float, f.readline().split()[:2])
                atoms_hdr = f.readline().strip()  # ITEM: ATOMS id type x y z vx vy vz ...
                cols = atoms_hdr.split()[2:]
                idx  = {c:i for i,c in enumerate(cols)}
                need = ["id","type","x","y","z","vx","vy","vz"]
                missing = [k for k in need if k not in idx]
                if missing:
                    raise SystemExit(f"Dump missing columns {missing} in {path}")
                ids   = np.empty(natoms, dtype=np.int64)
                types = np.empty(natoms, dtype=np.int32)
                pos   = np.empty((natoms,3), dtype=np.float64)
                vel   = np.empty((natoms,3), dtype=np.float64)
                for i in range(natoms):
                    parts   = f.readline().split()
                    ids[i]   = int(parts[idx["id"]])
                    types[i] = int(parts[idx["type"]])
                    pos[i,0] = float(parts[idx["x"]]);  pos[i,1] = float(parts[idx["y"]]);  pos[i,2] = float(parts[idx["z"]])
                    vel[i,0] = float(parts[idx["vx"]]); vel[i,1] = float(parts[idx["vy"]]); vel[i,2] = float(parts[idx["vz"]])
                box = (xlo,xhi,ylo,yhi,zlo,zhi)
                yield box, ids, types, pos, vel

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", required=True, help="runs/.../dumps/atoms.spike.*.lammpstrj")
    ap.add_argument("--r_core_nm", type=float, default=0.8, help="core cylinder radius [nm]")
    ap.add_argument("--frame_dt_ps", type=float, default=0.05, help="time spacing between frames [ps]")
    ap.add_argument("--out", required=True, help="output CSV path")
    args = ap.parse_args()

    r_core_A = args.r_core_nm * A_PER_NM

    Tl_list, t_list = [], []

    for f, (box, _ids, types, pos, vel) in enumerate(parse_sequence(args.pattern)):
        xlo,xhi,ylo,yhi,_,_ = box
        Lx, Ly = xhi-xlo, yhi-ylo
        cx, cy = xlo + 0.5*Lx, ylo + 0.5*Ly
        # radial mask
        r = np.sqrt((pos[:,0]-cx)**2 + (pos[:,1]-cy)**2)
        mask = r <= r_core_A
        n = int(mask.sum())
        t_list.append(f*args.frame_dt_ps)
        if n == 0:
            Tl_list.append(np.nan)
            continue
        m_amu = np.array([MASS_AMU[int(t)] for t in types[mask]], dtype=float)  # amu

        # Remove COM drift before computing kinetic temperature
        vcm   = (m_amu[:,None] * vel[mask]).sum(axis=0) / m_amu.sum()
        vrel  = vel[mask] - vcm
        v2    = np.einsum('ij,ij->i', vrel, vrel)               # (Å/ps)^2
        KE_eV = 0.5 * MVV2E * m_amu * v2                        # eV per atom
        KE_tot = KE_eV.sum()
        Tl = (2.0/3.0) * (KE_tot / (n * KB_eV_per_K))
        Tl_list.append(Tl)

    df = pd.DataFrame({"t_ps": np.asarray(t_list), "Tl_core_K": np.asarray(Tl_list)})
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, index=False)

    # Meta for reproducibility
    with open(args.out.replace(".csv",".meta.txt"), "w") as m:
        m.write(f"r_core_nm = {args.r_core_nm}\n")
        m.write(f"frame_dt_ps = {args.frame_dt_ps}\n")
        m.write("basis = atoms-only kinetic definition in metal units (COM-detrended)\n")

if __name__ == "__main__":
    main()
