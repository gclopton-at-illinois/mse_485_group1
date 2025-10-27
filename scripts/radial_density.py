#!/usr/bin/env python3
import argparse, os, glob, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MASS_AMU = {1: 140.116, 2: 15.999}
AMU_TO_G = 1.66053906660e-24
A3_TO_CM3 = 1e-24  # 1 Å^3 = 1e-24 cm^3

def read_frames(pattern):
    paths = sorted(glob.glob(pattern))
    for p in paths:
        with open(p, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if line.startswith("ITEM: TIMESTEP"):
                    # timestep line, then number
                    _ts = int(f.readline().strip())
                    _ = f.readline()             # ITEM: NUMBER OF ATOMS
                    natoms = int(f.readline())
                    _ = f.readline()             # ITEM: BOX BOUNDS ...
                    xlo,xhi = map(float, f.readline().split()[:2])
                    ylo,yhi = map(float, f.readline().split()[:2])
                    zlo,zhi = map(float, f.readline().split()[:2])
                    _ = f.readline()             # ITEM: ATOMS ...
                    cols = f.readline().split()
                    # The first data row already consumed! Fix by reading natoms-1 more, but we lost one row.
                    # Safer: re-open and parse with a small state machine per block.
                # If we hit a dump with multiple frames, better to parse block-wise.

def _step_key(path: str) -> int:
    m = re.search(r"(\d+)(?:\D*)$", os.path.basename(path))
    return int(m.group(1)) if m else 0

def parse_lammpstrj(pattern):
    """Return a list of (box, ids, types, pos, vel) per frame across all files matching pattern in numeric order."""
    frames = []
    for path in sorted(glob.glob(pattern), key=_step_key):
        with open(path, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if not line.startswith("ITEM: TIMESTEP"):
                    continue
                timestep = int(f.readline().strip())
                assert f.readline().startswith("ITEM: NUMBER OF ATOMS")
                natoms = int(f.readline().strip())
                bounds_hdr = f.readline().strip()
                assert bounds_hdr.startswith("ITEM: BOX BOUNDS")
                xlo,xhi = map(float, f.readline().split()[:2])
                ylo,yhi = map(float, f.readline().split()[:2])
                zlo,zhi = map(float, f.readline().split()[:2])
                atoms_hdr = f.readline().strip()
                assert atoms_hdr.startswith("ITEM: ATOMS")
                atom_cols = atoms_hdr.split()[2:]
                # build index map
                col_idx = {c:i for i,c in enumerate(atom_cols)}
                need = ["id","type","x","y","z"]
                has_v = all(k in col_idx for k in ["vx","vy","vz"])
                ids = np.empty(natoms, dtype=np.int64)
                types = np.empty(natoms, dtype=np.int32)
                pos = np.empty((natoms,3), dtype=np.float64)
                vel = np.empty((natoms,3), dtype=np.float64) if has_v else None
                for i in range(natoms):
                    parts = f.readline().split()
                    ids[i]   = int(parts[col_idx["id"]])
                    types[i] = int(parts[col_idx["type"]])
                    pos[i,0] = float(parts[col_idx["x"]]); pos[i,1] = float(parts[col_idx["y"]]); pos[i,2] = float(parts[col_idx["z"]])
                    if has_v:
                        vel[i,0] = float(parts[col_idx["vx"]]); vel[i,1] = float(parts[col_idx["vy"]]); vel[i,2] = float(parts[col_idx["vz"]])
                box = (xlo,xhi,ylo,yhi,zlo,zhi)
                frames.append((box,ids,types,pos,vel))
    return frames

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", required=True)
    ap.add_argument("--bin_A", type=float, default=1.0)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    frames = parse_lammpstrj(args.pattern)
    if len(frames) == 0:
        raise SystemExit(f"No frames matched pattern {args.pattern}")

    # Geometry from first frame
    (xlo,xhi,ylo,yhi,zlo,zhi), ids0, types0, pos0, _ = frames[0]
    Lx, Ly, Lz = xhi-xlo, yhi-ylo, zhi-zlo
    cx, cy = xlo+0.5*Lx, ylo+0.5*Ly
    r_max = 0.5*min(Lx, Ly) - args.bin_A
    nbins = max(1, int(np.floor(r_max/args.bin_A)))
    r_edges = np.linspace(0.0, nbins*args.bin_A, nbins+1)
    r_centers = 0.5*(r_edges[:-1] + r_edges[1:])
    shell_vol = 2.0*np.pi*r_centers*args.bin_A*Lz     # Å^3 per shell

    def rho_frame(box, types, pos):
        r = np.sqrt((pos[:,0]-cx)**2 + (pos[:,1]-cy)**2)
        masses_g = np.vectorize(MASS_AMU.get)(types) * AMU_TO_G
        mass_sum_g, _ = np.histogram(r, bins=r_edges, weights=masses_g)
        return mass_sum_g / (shell_vol * A3_TO_CM3)   # g/cm^3

    RHO = np.vstack([rho_frame(*fr[:1], fr[2], fr[3]) for fr in frames])  # (T,R)
    far = max(1, int(0.9*nbins))
    rho0 = RHO[0, far:].mean()
    threshold = 0.9*rho0

    R_track = []
    for row in RHO:
        mask = row < threshold
        if mask.any():
            idx = int(np.argmax(mask))
            R_track.append(r_centers[idx])
        else:
            R_track.append(np.nan)
    R_track = np.array(R_track)

    # Save CSVs
    pd.DataFrame({"frame": np.arange(len(frames)), "R_track_A": R_track}).to_csv(
        os.path.join(args.outdir, "R_track_vs_frame.csv"), index=False
    )
    df = pd.DataFrame(RHO, columns=[f"{rc:.3f}" for rc in r_centers])
    df.insert(0, "frame", np.arange(len(frames)))
    df.to_csv(os.path.join(args.outdir, "rho_r_vs_frame.csv"), index=False)
    with open(os.path.join(args.outdir, "rho_meta.txt"), "w") as fh:
        fh.write(f"rho0_farfield_g_per_cm3 = {rho0:.6f}\n")
        fh.write(f"threshold = {threshold:.6f}\n")
        fh.write(f"r_bin_A = {args.bin_A}\n")
        fh.write(f"Lz_A = {Lz}\n")
        fh.write(f"nbins = {nbins}\n")

    # Quick-look figure
    fig, ax = plt.subplots()
    T = len(frames)
    picks = sorted(set([0, max(1,T//3), max(2,2*T//3), T-1]))
    for f in picks:
        ax.plot(r_centers/10.0, RHO[f], label=f"frame {f}")
        if not np.isnan(R_track[f]):
            ax.axvline(R_track[f]/10.0, ls="--", alpha=0.5)
    ax.axhline(threshold, ls=":", alpha=0.7, label="0.9·rho0")
    ax.set_xlabel("r [nm]"); ax.set_ylabel("mass density [g/cm³]"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "rho_profiles.png"), dpi=180)

if __name__ == "__main__":
    main()
