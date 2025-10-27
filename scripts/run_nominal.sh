#!/usr/bin/env bash
# Phase 0 nominal shot runner — CeO2 TTM–MD
# Usage:
#   scripts/run_nominal.sh [--lmp /path/to/lmp] [--run-dir runs/phase0_nominal/<stamp>] [--dry-run]
# Requirements:
#   - Python with PyYAML in your conda env
#   - Repo layout as documented (inputs/, config/, etc.)

set -Eeuo pipefail
trap 'echo "[ERR] line $LINENO: $BASH_COMMAND" >&2' ERR

# --------------- CLI ---------------
LMP_BIN="${LMP:-lmp}"         # override with env LMP or --lmp
RUN_DIR=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --lmp)       LMP_BIN="$2"; shift 2 ;;
    --run-dir)   RUN_DIR="$2"; shift 2 ;;
    --dry-run)   DRY_RUN=1; shift ;;
    -h|--help)
      echo "Usage: $0 [--lmp /path/to/lmp] [--run-dir runs/phase0_nominal/<stamp>] [--dry-run]"
      exit 0
      ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# --------------- Parse physics.yaml (and pick up TTM paths) ---------------
IFS= read -r -d '' PY_HELPER << 'PYCODE' || :
import sys, yaml, pathlib
p = pathlib.Path("config/physics.yaml")
d = yaml.safe_load(p.read_text())

def req(dic, *path):
    x = dic
    for k in path:
        if k not in x:
            sys.stderr.write(f"Missing key in config/physics.yaml: {'.'.join(path)}\n")
            sys.exit(3)
        x = x[k]
    return x

units_family   = req(d,'units','family')
dt_fs          = float(req(d,'time','dt_fs'))
t_total_ps     = float(req(d,'time','t_total_ps'))
bx             = req(d,'boundaries','x')
by             = req(d,'boundaries','y')
bz             = req(d,'boundaries','z')

r0_nm          = float(req(d,'source','r0_nm'))
Se_eV_per_A    = float(req(d,'source','Se_eV_per_A'))

g0             = float(req(d,'ttm','coupling','g0_W_m3K'))
ke_file        = req(d,'ttm','ke_curve_file')
Ce_file        = req(d,'ttm','Ce_curve_file')
grid_file      = req(d,'ttm','grid_file')
src_file       = req(d,'source','pulse') and 'inputs/ttm/source_profile.yaml'

rim = d.get('forcefield',{}).get('rim_langevin',{})
gamma_ps_inv   = rim.get('gamma_ps_inv', None)
rim_width_nm   = rim.get('width_nm', None)
rim_target_K   = rim.get('target_K', None)

data_path      = req(d,'structure','data_path')

dt_ps = dt_fs * 1e-3
nsteps = int(round(t_total_ps / dt_ps))

def shquote(s):
    return "'" + str(s).replace("'", "'\"'\"'") + "'"

assign = {
  "UNITS": units_family,
  "DT": f"{dt_ps:.12g}",
  "NSTEPS": str(nsteps),
  "R0": f"{r0_nm:.12g}",
  "SE": f"{Se_eV_per_A:.12g}",
  "G0": f"{g0:.12g}",
  "BX": bx, "BY": by, "BZ": bz,
  "DATAFILE": data_path,
  "KE_FILE": ke_file,
  "CE_FILE": Ce_file,
  "GRID_FILE": grid_file,
  "SRC_FILE": src_file,
}

if gamma_ps_inv is not None: assign["RIM_GAMMA"] = f"{float(gamma_ps_inv):.12g}"
if rim_width_nm  is not None: assign["RIM_WIDTH_NM"] = f"{float(rim_width_nm):.12g}"
if rim_target_K  is not None: assign["RIM_TK"] = f"{float(rim_target_K):.12g}"

for k,v in assign.items():
    print(f"{k}={shquote(v)}")
PYCODE

# shellcheck disable=SC2046
eval "$(
  python3 - <<PYCODE
$PY_HELPER
PYCODE
)"

# --------------- Resolve run directory & snapshot ---------------
timestamp="$(date +%Y%m%d_%H%M%S)"
if [[ -z "$RUN_DIR" ]]; then
  RUN_DIR="runs/phase0_nominal/${timestamp}"
fi

SNAP="$RUN_DIR/inputs_snapshot"
mkdir -p "$RUN_DIR"/{logs,dumps,post} \
         "$SNAP"/{lammps/includes,structure,ttm} \
         "$RUN_DIR"

# Copy snapshot of configs and inputs for provenance
cp -f config/physics.yaml              "$SNAP/physics.yaml"        || true
cp -f config/metrics.yaml              "$SNAP/metrics.yaml"        || true
cp -f config/figure_style.yaml         "$SNAP/figure_style.yaml"   || true

cp -f "$DATAFILE"                      "$SNAP/structure/"
cp -f inputs/lammps/in.main.lmp        "$SNAP/lammps/"
cp -f inputs/lammps/in.equil.lmp       "$SNAP/lammps/"
cp -f inputs/lammps/in.ttm_spike.lmp   "$SNAP/lammps/"
cp -f inputs/lammps/includes/*.in      "$SNAP/lammps/includes/"    2>/dev/null || true
[[ -f inputs/lammps/tables/zbl.table ]] && cp -f inputs/lammps/tables/zbl.table "$SNAP/lammps/"

cp -f "$KE_FILE"                       "$SNAP/ttm/" || true
cp -f "$CE_FILE"                       "$SNAP/ttm/" || true
cp -f "$GRID_FILE"                     "$SNAP/ttm/" || true
cp -f "$SRC_FILE"                      "$SNAP/ttm/" || true

# --------------- Minimal manifest ---------------
cat > "$RUN_DIR/manifest.yaml" <<EOF
phase: 0
units: ${UNITS}
dt_ps: ${DT}
nsteps: ${NSTEPS}
source:
  r0_nm: ${R0}
  Se_eV_per_A: ${SE}
ttm:
  g0_W_m3K: ${G0}
paths:
  datafile: ${DATAFILE}
  ke_curve_file: ${KE_FILE}
  Ce_curve_file: ${CE_FILE}
  grid_file: ${GRID_FILE}
  src_file: ${SRC_FILE}
boundaries: { x: ${BX}, y: ${BY}, z: ${BZ} }
snapshot_dir: ${SNAP}
timestamp: ${timestamp}
EOF

# ---------- Precompute Te.in (cylindrical Gaussian on TTM grid) ----------
TEF="${RUN_DIR}/post/Te.in"
python3 - <<PY || { echo "ERROR: Te.in generation failed" >&2; exit 88; }
import os, math, sys

DATAFILE = r"""${DATAFILE}"""
R0_nm    = float(r"""${R0}""")
TE_BASE  = 300.0
TE_CORE  = 20000.0
grid_A   = 1.0  # Å target spacing (Phase 0)

# Parse xlo/xhi etc. from LAMMPS data file
xlo=xhi=ylo=yhi=zlo=zhi=None
with open(DATAFILE,"r") as f:
    for ln in f:
        s=ln.split()
        if len(s)==4 and s[2:]==["xlo","xhi"]:
            xlo,xhi = float(s[0]), float(s[1])
        elif len(s)==4 and s[2:]==["ylo","yhi"]:
            ylo,yhi = float(s[0]), float(s[1])
        elif len(s)==4 and s[2:]==["zlo","zhi"]:
            zlo,zhi = float(s[0]), float(s[1])
        if None not in (xlo,xhi,ylo,yhi,zlo,zhi):
            break
if None in (xlo,xhi,ylo,yhi,zlo,zhi):
    print("Could not parse box from data file", file=sys.stderr); sys.exit(2)

Lx, Ly, Lz = xhi-xlo, yhi-ylo, zhi-zlo
NX = int(math.ceil(Lx/grid_A));  NY = int(math.ceil(Ly/grid_A));  NZ = int(math.ceil(Lz/grid_A))
dx, dy = Lx/NX, Ly/NY
sigma = (R0_nm*10.0)/math.sqrt(2.0)  # Å

tef = r"""${TEF}"""
os.makedirs(os.path.dirname(tef), exist_ok=True)
with open(tef, "w") as out:
    for iz in range(1, NZ+1):
        for iy in range(1, NY+1):
            y = (iy-0.5)*dy - 0.5*Ly
            for ix in range(1, NX+1):
                x = (ix-0.5)*dx - 0.5*Lx
                r2 = x*x + y*y
                Te = TE_BASE + (TE_CORE-TE_BASE)*math.exp(-r2/(2.0*sigma*sigma))
                out.write(f"{ix} {iy} {iz} {Te:.6f}\n")

# tiny manifest
with open(tef + ".meta.txt","w") as m:
    m.write(f"NX NY NZ = {NX} {NY} {NZ}\n")
    m.write(f"expected_lines = {NX*NY*NZ}\n")
    m.write(f"Lx Ly Lz (A) = {Lx} {Ly} {Lz}\n")
    m.write(f"R0_nm = {R0_nm}\n")
PY

# Quick sanity print to job logs
wc -l "$TEF" || true
head -n 2 "$TEF" || true
tail -n 2 "$TEF" || true

# --------------- Build LAMMPS command ---------------
LAMMPS_IN=inputs/lammps/in.main.lmp
LOGFILE="${RUN_DIR}/logs/log.lammps"
SCREEN="${RUN_DIR}/logs/screen.out"

mkdir -p "$(dirname "$SCREEN")" "$(dirname "$LOGFILE")"
: > "$SCREEN"
: > "$LOGFILE"

LAMMPS_CMD=(
  "$LMP_BIN" -echo both
  -in "$LAMMPS_IN" -log "$LOGFILE"
  -var UNITS "$UNITS" -var DATAFILE "$DATAFILE"
  -var DT "$DT" -var NSTEPS "$NSTEPS"
  -var R0 "$R0" -var SE "$SE" -var G0 "$G0"
  -var KE_FILE "$KE_FILE" -var CE_FILE "$CE_FILE"
  -var GRID_FILE "$GRID_FILE" -var SRC_FILE "$SRC_FILE"
  -var BX "$BX" -var BY "$BY" -var BZ "$BZ"
  -var VARS_FILE "$RUN_DIR/vars.lmp"
  -var OUTDIR "$RUN_DIR"
  -var TEF "$TEF"
)

# Optional rim parameters (only if present in YAML and consumed by your includes)
[[ ${RIM_GAMMA:-}    ]] && LAMMPS_CMD+=( -var RIM_GAMMA "$RIM_GAMMA" )
[[ ${RIM_WIDTH_NM:-} ]] && LAMMPS_CMD+=( -var RIM_WIDTH_NM "$RIM_WIDTH_NM" )
[[ ${RIM_TK:-}       ]] && LAMMPS_CMD+=( -var RIM_TK "$RIM_TK" )

echo "=== Phase 0 nominal shot ==="
echo "Run dir:    $RUN_DIR"
echo "LMP binary: $LMP_BIN"
echo "Input deck: $LAMMPS_IN"
echo "DT (ps):    $DT"
echo "NSTEPS:     $NSTEPS"
echo "R0 (nm):    $R0"
echo "Se (eV/Å):  $SE"
echo "g0 (W/m3K): $G0"
echo

if [[ $DRY_RUN -eq 1 ]]; then
  printf '%q ' "${LAMMPS_CMD[@]}"; echo
  echo "(dry-run: not executing)"
  exit 0
fi

# --------------- Launch ---------------
LAUNCH=()
if command -v srun >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
  LAUNCH=(srun -u -n 1)
fi

STD_BUF=()
if command -v stdbuf >/dev/null 2>&1; then
  STD_BUF=(stdbuf -oL -eL)
fi

"${LAUNCH[@]}" "${STD_BUF[@]}" "${LAMMPS_CMD[@]}" >"$SCREEN" 2>&1
RC=$?

if (( RC != 0 )); then
  echo "LAMMPS exited with RC=$RC" >&2
  echo "---- screen.out (head) ----" >&2
  sed -n '1,80p'   "$SCREEN" >&2 || true
  echo "---- screen.out (tail) ----" >&2
  tail -n 80       "$SCREEN" >&2 || true
  echo "---- log.lammps (tail) ----" >&2
  tail -n 120      "$LOGFILE" >&2 || true
  exit "$RC"
fi

echo "Done. Logs at: $RUN_DIR/logs"
echo "Snapshot at:   $SNAP"
