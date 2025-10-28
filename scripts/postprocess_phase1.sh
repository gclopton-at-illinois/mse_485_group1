#!/usr/bin/env bash
# Phase 1 postprocessing: rho(r,t)+R_track(t), WS-defects vs time, Te/Tl overlay, energy balance.
set -Eeuo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ---------------- Resolve RUN dir ----------------
# Prefer explicit RUN=..., else runs/phase1_baseline/LATEST, else newest under common roots
if [[ -n "${RUN:-}" && -d "${RUN:-/dev/null}" ]]; then
  : # keep as provided
else
  CAND=""
  if [[ -L "$REPO_ROOT/runs/phase1_baseline/LATEST" ]]; then
    CAND=$(readlink -f "$REPO_ROOT/runs/phase1_baseline/LATEST" || true)
    [[ -d "$CAND" ]] || CAND=""
  fi
  if [[ -z "$CAND" ]]; then
    CAND=$(ls -1dt \
      "$REPO_ROOT/runs/phase1_baseline"/*/ \
      /u/gclopton/scratch/lammps/runs/phase1_baseline/*/ \
      /scratch/gclopton/lammps/runs/phase1_baseline/*/ 2>/dev/null | head -1 || true)
  fi
  RUN="$CAND"
fi

[[ -n "${RUN:-}" && -d "${RUN:-/dev/null}" ]] || { echo "ERROR: cannot resolve RUN (try: RUN=/path/to/run scripts/postprocess_phase1.sh)"; exit 1; }

echo "[p1-post] RUN=${RUN}"

DUMPS="$RUN/dumps"
POST="$RUN/post"
FIGS="$RUN/figures"
mkdir -p "$POST" "$FIGS"

# ---------------- Frame timing (fallback if unknown) ----------------
# If you know the actual stride, export FRAME_DT_PS before calling; else default 0.05 ps
FRAME_DT_PS=${FRAME_DT_PS:-0.05}

# ---------------- 1) Radial density & R_track(t) ----------------
# Outputs (by convention):
#   $POST/rho_r_vs_frame.csv
#   $POST/R_track_vs_frame.csv
#   $POST/rho_meta.txt  (may include threshold, bin size, etc.)
python "$REPO_ROOT/scripts/radial_density.py" \
  --pattern "$DUMPS/atoms.spike.*.lammpstrj" \
  --bin_A 2.0 \
  --outdir "$POST"

# ---------------- 2) Defect density vs time (Wigner–Seitz proxy) ----------------
FIRST=$(ls "$DUMPS"/atoms.spike.*.lammpstrj 2>/dev/null | sort -V | head -1 || true)
[[ -n "$FIRST" && -f "$FIRST" ]] || { echo "ERROR: no spike frames found under $DUMPS"; exit 2; }

python "$REPO_ROOT/scripts/ws_defects.py" \
  --pattern "$DUMPS/atoms.spike.*.lammpstrj" \
  --ref "$FIRST" \
  --frame_dt_ps "$FRAME_DT_PS" \
  --disp_thr_A 0.6 \
  --cyl_radius_nm 8.0 \
  --out "$POST/ndef_vs_time.csv"

# ---------------- 3) Core Tl (and optional Te if script supports it) ----------------
# tl_core.py currently writes Tl_core.csv. If you later extend it to also write Te_core.csv,
# make_phase1_figs.py will automatically overlay both.
python "$REPO_ROOT/scripts/tl_core.py" \
  --pattern "$DUMPS/atoms.spike.*.lammpstrj" \
  --r_core_nm 0.8 \
  --frame_dt_ps "$FRAME_DT_PS" \
  --out "$POST/Tl_core.csv"

# ---------------- 4) Figures (Phase-1 trio) ----------------
python "$REPO_ROOT/scripts/make_phase1_figs.py" \
  --post "$POST" \
  --figs "$FIGS"

# ---------------- 5) Energy accounting sanity gate ----------------
if [[ -f "$REPO_ROOT/scripts/energy_balance.py" ]]; then
  python "$REPO_ROOT/scripts/energy_balance.py" \
    --run "$RUN" \
    --physics "$REPO_ROOT/config/physics.yaml" \
    --tol 0.02 || true
fi

echo "[p1-post] Done. CSVs in $POST, figures in $FIGS"
