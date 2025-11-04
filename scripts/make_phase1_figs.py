#!/usr/bin/env python3

import argparse, os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def _safe_read_csv(path):
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)

def fig_rho_radial(post, figs):
    rho_path = os.path.join(post, "rho_r_vs_frame.csv")
    R_path   = os.path.join(post, "R_track_vs_frame.csv")
    rho = _safe_read_csv(rho_path)
    R   = _safe_read_csv(R_path)
    if rho is None or R is None:
        return  # nothing to do

    frames = rho["frame"].to_numpy()
    rA = np.array([float(c) for c in rho.columns[1:]], dtype=float)
    rnm = rA / 10.0

    # Choose representative frames: pre-pulse, ~1/3, ~2/3, final
    idxs = sorted(set([
        0,
        max(1, len(frames)//3),
        max(2, (2*len(frames))//3),
        len(frames)-1
    ]))

    # Threshold line from meta 
    thr = None
    meta = os.path.join(post, "rho_meta.txt")
    if os.path.exists(meta):
        with open(meta) as fh:
            for ln in fh:
                if "threshold" in ln and "=" in ln:
                    try:
                        thr = float(ln.split("=",1)[1].strip())
                    except Exception:
                        pass

    plt.figure()
    for f in idxs:
        y = rho.iloc[f,1:].to_numpy(dtype=float)
        label = f"frame {int(frames[f])}"
        plt.plot(rnm, y, lw=2, label=label)
        RtA = R["R_track_A"].iloc[f] if f < len(R) else np.nan
        if np.isfinite(RtA):
            plt.scatter([RtA/10.0],[np.interp(RtA, rA, y)], s=30, marker="o")
            plt.axvline(RtA/10.0, ls="--", alpha=0.4)

    if thr is not None:
        plt.axhline(thr, ls=":", alpha=0.7, label="0.9 · ρ₀")

    plt.xlabel("r [nm]")
    plt.ylabel("mass density [g/cm³]")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figs, "rho_radial_overlays.png"), dpi=200)

def fig_ndef_time(post, figs):
    path = os.path.join(post, "ndef_vs_time.csv")
    df = _safe_read_csv(path)
    if df is None:
        return
    # assume first numeric column after t_ps is the density
    t = df["t_ps"].to_numpy()
    y = df[df.columns[1]].to_numpy(dtype=float)

    peak = float(np.nanmax(y)) if y.size else np.nan
    # Optional metadata (t_half, t_1e)
    meta = os.path.splitext(path)[0] + ".meta.txt"
    t_half = t_1e = None
    if os.path.exists(meta):
        with open(meta) as fh:
            for ln in fh:
                if ln.startswith("t_half_ps"): t_half = float(ln.split("=")[1])
                if ln.startswith("t_1e_ps"):   t_1e   = float(ln.split("=")[1])

    # plateau (last 10% or ≥3 points)
    k = max(3, int(0.1*len(y)))
    plateau = float(np.nanmean(y[-k:])) if y.size else np.nan

    def cross_time(frac: float) -> float:
        if not np.isfinite(plateau): return np.nan
        thr = frac * plateau
        idx = np.where(y >= thr)[0]
        if idx.size == 0: return np.nan
        j = int(idx[0])
        if j == 0: return float(t[0])
        t0, t1 = t[j-1], t[j]
        y0, y1 = y[j-1], y[j]
        if y1 == y0: return float(t1)
        return float(t0 + (thr - y0) * (t1 - t0) / (y1 - y0))

    t50, t90 = cross_time(0.5), cross_time(0.9)

    plt.figure()
    plt.plot(t, y, lw=2)
    if np.isfinite(peak):
        plt.axhline(peak, ls=":", alpha=0.5, label="peak")
    if (t_half is not None) and np.isfinite(t_half):
        plt.axvline(t_half, ls="--", alpha=0.7, label=r"$t_{1/2}$")
    if (t_1e is not None) and np.isfinite(t_1e):
        plt.axvline(t_1e,  ls="--", alpha=0.7, label=r"$t_{1/e}$")
    if np.isfinite(t50):
        plt.axvline(t50, ls="--", alpha=0.8, label=f"t50={t50:.2f} ps")
    if np.isfinite(t90):
        plt.axvline(t90, ls="--", alpha=0.8, label=f"t90={t90:.2f} ps")

    plt.xlabel("t [ps]")
    plt.ylabel(r"defect density [nm$^{-3}$]")
    handles, labels = plt.gca().get_legend_handles_labels()
    if labels:
        plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(os.path.join(figs, "defects_time_series.png"), dpi=200)

def fig_TeTl_time(post, figs):
    tl = _safe_read_csv(os.path.join(post, "Tl_core.csv"))
    te = _safe_read_csv(os.path.join(post, "Te_core.csv"))  # optional

    if tl is None and te is None:
        return

    plt.figure()
    if tl is not None:
        plt.plot(tl["t_ps"].to_numpy(), tl["Tl_core_K"].to_numpy(), lw=2, label=r"$T_l$ (core)")
        t_ref = tl["t_ps"].to_numpy()
    else:
        t_ref = None

    if te is not None:
        plt.plot(te["t_ps"].to_numpy(), te["Te_core_K"].to_numpy(), lw=2, label=r"$T_e$ (core)")

    plt.xlabel("t [ps]")
    plt.ylabel("Temperature [K]")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figs, "Te_Tl_core.png"), dpi=200)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--post", required=True, help="runs/.../post directory")
    ap.add_argument("--figs", required=True, help="runs/.../figures directory")
    args = ap.parse_args()
    os.makedirs(args.figs, exist_ok=True)
    fig_rho_radial(args.post, args.figs)
    fig_ndef_time(args.post, args.figs)
    fig_TeTl_time(args.post, args.figs)

if __name__ == "__main__":
    main()
