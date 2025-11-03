#!/usr/bin/env python3
"""
Radial mass density ρ(r,t) in cylindrical shells around the box center,
+ R_track(t) pickoff as the first r where ρ < 0.9*ρ0 (far-field baseline).

Inputs
------
--pattern runs/.../dumps/atoms.spike.*.lammpstrj
  Expect frames with: id type x y z [vx vy vz] (vels not required here)

Options
-------
--bin_A <float>                Radial bin width in Å (default 2.0)
--smooth_window_bins <int>     Moving-average window over r (odd, default 5)
--drop_inner_bins <int>        Skip this many innermost bins in the *plot only* (CSV unaffected; default 0)
--pickoff_rmin_nm <float>      Ignore r < this when picking first ρ < 0.9·ρ0 (default 0.0 nm)
--outdir <path>                Output directory (CSV + quick-look PNG)

Outputs
-------
- rho_r_vs_frame.csv           rows: frames; cols: 'frame', r_centers [Å] as headers
- R_track_vs_frame.csv         cols: frame, R_track_A
- rho_meta.txt                 baseline ρ0, threshold, Lz, bins, smoothing
- rho_profiles.png             quick-look overlays with pickoff markers
"""
import argparse, os, glob, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MASS_AMU = {1: 140.116, 2: 15.999}
AMU_TO_G = 1.66053906660e-24
A3_TO_CM3 = 1e-24  # 1 Å^3 = 1e-24 cm^3

def _step_key(path: str) -> int:
    m = re.search(r"(\d+)(?:\D*)$", os.path.basename(path))
    return int(m.group(1)) if m else 0

def parse_lammpstrj(pattern):
    frames = []
    for path in sorted(glob.glob(pattern), key=_step_key):
        with open(path, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if not line.startswith("ITEM: TIMESTEP"):
                    continue
                _ = int(f.readline().strip())                       # timestep
                assert f.readline().startswith("ITEM: NUMBER OF ATOMS")
                natoms = int(f.readline().strip())
                bounds_hdr = f.readline().strip()
                assert bounds_hdr.startswith("ITEM: BOX BOUNDS")
                xlo,xhi = map(float, f.readline().split()[:2])
                ylo,yhi = map(float, f.readline().split()[:2])
                zlo,zhi = map(float, f.readline().split()[:2])
                atoms_hdr = f.readline().strip()
                assert atoms_hdr.startswith("ITEM: ATOMS")
                cols = atoms_hdr.split()[2:]
                idx = {c:i for i,c in enumerate(cols)}
                need = ["id","type","x","y","z"]
                missing = [k for k in need if k not in idx]
                if missing:
                    raise SystemExit(f"Dump missing columns {missing} in {path}")
                ids   = np.empty(natoms, dtype=np.int64)
                types = np.empty(natoms, dtype=np.int32)
                pos   = np.empty((natoms,3), dtype=np.float64)
                for i in range(natoms):
                    parts   = f.readline().split()
                    ids[i]   = int(parts[idx["id"]])
                    types[i] = int(parts[idx["type"]])
                    pos[i,0] = float(parts[idx["x"]]); pos[i,1] = float(parts[idx["y"]]); pos[i,2] = float(parts[idx["z"]])
                frames.append(((xlo,xhi,ylo,yhi,zlo,zhi), ids, types, pos))
    return frames

def moving_average(y: np.ndarray, w: int) -> np.ndarray:
    if w <= 1 or w % 2 == 0:  # enforce odd window; 1 disables smoothing
        return y
    pad = w // 2
    ypad = np.pad(y, (pad, pad), mode="reflect")
    kern = np.ones(w, dtype=float) / w
    return np.convolve(ypad, kern, mode="valid")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", required=True)
    ap.add_argument("--bin_A", type=float, default=2.0)
    ap.add_argument("--smooth_window_bins", type=int, default=5,
                    help="odd window length; 1 disables smoothing")
    ap.add_argument("--drop_inner_bins", type=int, default=0,
                    help="skip this many inner bins in the quick-look plot only")
    ap.add_argument("--pickoff_rmin_nm", type=float, default=0.0,
                    help="ignore r < this value when searching for ρ < 0.9·ρ0")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    frames = parse_lammpstrj(args.pattern)
    if len(frames) == 0:
        raise SystemExit(f"No frames matched pattern {args.pattern}")

    # Geometry from first frame
    (xlo,xhi,ylo,yhi,zlo,zhi), ids0, types0, pos0 = frames[0]
    Lx, Ly, Lz = xhi-xlo, yhi-ylo, zhi-zlo
    cx, cy = xlo+0.5*Lx, ylo+0.5*Ly
    r_max = 0.5*min(Lx, Ly) - args.bin_A
    nbins = max(1, int(np.floor(r_max/args.bin_A)))
    r_edges = np.linspace(0.0, nbins*args.bin_A, nbins+1)
    r_centers = 0.5*(r_edges[:-1] + r_edges[1:])
    r_in  = r_edges[:-1]
    r_out = r_edges[1:]
    shell_vol = np.pi * (r_out**2 - r_in**2) * Lz

    def rho_frame(types, pos):
        r = np.sqrt((pos[:,0]-cx)**2 + (pos[:,1]-cy)**2)
        masses_g = np.array([MASS_AMU[int(t)] for t in types], dtype=float) * AMU_TO_G
        mass_sum_g, _ = np.histogram(r, bins=r_edges, weights=masses_g)
        rho = mass_sum_g / (shell_vol * A3_TO_CM3)   # g/cm^3
        return rho

    RHO = np.vstack([rho_frame(fr[2], fr[3]) for fr in frames])  # (T,R)

    # Light, time-independent smoothing over r (keeps elbow location stable)
    if args.smooth_window_bins > 1 and args.smooth_window_bins % 2 == 1:
        RHO = np.vstack([moving_average(row, args.smooth_window_bins) for row in RHO])

    # Far-field baseline ρ0: average of last max(5, 10%) bins at t=0
    tail = max(5, max(1, nbins//10))
    rho0 = float(np.nanmean(RHO[0, -tail:]))
    threshold = 0.9 * rho0

    # R_track pickoff (optionally ignore inner-most region)
    start_idx = int(np.searchsorted(r_centers/10.0, args.pickoff_rmin_nm))
    R_track = []
    for row in RHO:
        below = np.nonzero(row[start_idx:] < threshold)[0]
        R_track.append(r_centers[start_idx + below[0]] if below.size else np.nan)
    R_track = np.asarray(R_track, dtype=float)

    # Save CSVs
    pd.DataFrame({"frame": np.arange(len(frames)), "R_track_A": R_track}).to_csv(
        os.path.join(args.outdir, "R_track_vs_frame.csv"), index=False
    )
    df = pd.DataFrame(RHO, columns=[f"{rc:.3f}" for rc in r_centers])
    df.insert(0, "frame", np.arange(len(frames)))
    df.to_csv(os.path.join(args.outdir, "rho_r_vs_frame.csv"), index=False)

    # Meta
    with open(os.path.join(args.outdir, "rho_meta.txt"), "w") as fh:
        fh.write(f"rho0_farfield_g_per_cm3 = {rho0:.6f}\n")
        fh.write(f"threshold = {threshold:.6f}\n")
        fh.write(f"r_bin_A = {args.bin_A}\n")
        fh.write(f"smooth_window_bins = {args.smooth_window_bins}\n")
        fh.write(f"drop_inner_bins = {args.drop_inner_bins}\n")
        fh.write(f"pickoff_rmin_nm = {args.pickoff_rmin_nm}\n")
        fh.write(f"Lz_A = {Lz}\n")
        fh.write(f"nbins = {nbins}\n")

    # Quick-look figure (optionally drop inner bins from the plot only)
    fig, ax = plt.subplots()
    T = len(frames)
    picks = sorted(set([0, max(1, T//3), max(2, 2*T//3), T-1]))
    plot_start = max(0, int(args.drop_inner_bins))

    for f in picks:
        ax.plot((r_centers/10.0)[plot_start:], RHO[f, plot_start:], lw=2, label=f"frame {f}")
        if not np.isnan(R_track[f]):
            # Only draw marker/line if the pickoff radius is in the plotted range
            if plot_start < len(r_centers) and R_track[f] >= r_centers[plot_start]:
                j = int(np.argmin(np.abs(r_centers - R_track[f])))
                ax.scatter([R_track[f]/10.0], [RHO[f, j]], s=28)
                ax.axvline(R_track[f]/10.0, ls="--", alpha=0.5)

    ax.axhline(threshold, ls=":", alpha=0.7, label="0.9·ρ0")
    ax.set_xlabel("r [nm]")
    ax.set_ylabel("mass density [g/cm³]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "rho_profiles.png"), dpi=200)

if __name__ == "__main__":
    main()
