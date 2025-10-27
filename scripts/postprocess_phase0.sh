#!/usr/bin/env bash
set -Eeuo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Prefer explicit RUN=..., else try LATEST, else newest under /u or /scratch
if [[ -n "${RUN:-}" && -d "${RUN:-/dev/null}" ]]; then
  : # keep RUN as provided
else
  CAND=""
  if [[ -L "$REPO_ROOT/runs/phase0_nominal/LATEST" ]]; then
    CAND=$(readlink -f "$REPO_ROOT/runs/phase0_nominal/LATEST" || true)
    [[ -d "$CAND" ]] || CAND=""
  fi
  if [[ -z "$CAND" ]]; then
    CAND=$(ls -1dt /u/gclopton/scratch/lammps/runs/phase0_nominal/*/ /scratch/gclopton/lammps/runs/phase0_nominal/*/ 2>/dev/null | head -1 || true)
  fi
  RUN="$CAND"
fi

[[ -n "${RUN:-}" && -d "${RUN:-/dev/null}" ]] || { echo "ERROR: cannot resolve RUN (try: RUN=/path/to/run scripts/postprocess_phase0.sh)"; exit 1; }

echo "Postprocessing RUN=${RUN}"

DUMPS="$RUN/dumps"
POST="$RUN/post"
FIGS="$RUN/figures"
mkdir -p "$POST" "$FIGS"

# Spike dump cadence: 100 steps at DT=0.0005 ps -> 0.05 ps between frames
FRAME_DT_PS=${FRAME_DT_PS:-0.05}

# 1) Radial density & R_track(t)  (pure NumPy)  — wider bins to reduce noise
python "$REPO_ROOT/scripts/radial_density.py" \
  --pattern "$DUMPS/atoms.spike.*.lammpstrj" \
  --bin_A 2.0 \
  --outdir "$POST"

# 2) Defect density via displacement proxy (Phase-0: reference = first spike frame)
#    This makes n_def(t0)=0 and avoids the artificial early dip.
FIRST=$(ls "$DUMPS"/atoms.spike.*.lammpstrj 2>/dev/null | sort -V | head -1 || true)
[[ -n "$FIRST" && -f "$FIRST" ]] || { echo "ERROR: no spike frames found under $DUMPS"; exit 2; }

python "$REPO_ROOT/scripts/ws_defects.py" \
  --pattern "$DUMPS/atoms.spike.*.lammpstrj" \
  --ref     "$FIRST" \
  --frame_dt_ps "$FRAME_DT_PS" \
  --disp_thr_A 0.6 \
  --out     "$POST/ndef_vs_time.csv"

# 3) Core Tl(t) — pure NumPy (velocities)
python "$REPO_ROOT/scripts/tl_core.py" \
  --pattern "$DUMPS/atoms.spike.*.lammpstrj" \
  --r_core_nm 0.8 \
  --frame_dt_ps "$FRAME_DT_PS" \
  --out "$POST/Tl_core.csv"

# 4) Make the three Phase-0 figures
python "$REPO_ROOT/scripts/make_phase0_figs.py" \
  --post "$POST" \
  --figs "$FIGS"

echo "Done. CSVs in $POST, figures in $FIGS"
