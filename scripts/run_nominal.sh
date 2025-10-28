#!/usr/bin/env bash
# Phase 1 baseline runner — CeO2 TTM–MD (δ=0, g=g0)
# Usage:
#   scripts/run_nominal.sh [--lmp /path/to/lmp] [--run-dir runs/phase1_baseline/<stamp>] [--config config/physics.yaml] [--dry-run]

set -Eeuo pipefail
trap 'echo "[ERR] line $LINENO: $BASH_COMMAND" >&2' ERR

# -------- CLI --------
LMP_BIN="${LMP:-lmp}"
RUN_DIR=""
DRY_RUN=0
PHYS_YAML="config/physics.yaml"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --lmp)       LMP_BIN="$2"; shift 2 ;;
    --run-dir)   RUN_DIR="$2"; shift 2 ;;
    --config|-c) PHYS_YAML="$2"; shift 2 ;;
    --dry-run)   DRY_RUN=1;    shift ;;
    -h|--help)
      echo "Usage: $0 [--lmp PATH] [--run-dir DIR] [--config config/physics.yaml] [--dry-run]"
      exit 0
      ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

[[ -x "$LMP_BIN" ]] || { echo "ERROR: lmp not executable at $LMP_BIN" >&2; exit 127; }

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
cd "$REPO_ROOT"

# -------- Parse physics.yaml and emit shell assignments --------
eval "$(
  python3 - "$PHYS_YAML" <<'PYCODE'
import sys, yaml, pathlib, time, random
p = pathlib.Path(sys.argv[1])
d = yaml.safe_load(p.read_text())

def req(dic, *keys):
    cur = dic
    for k in keys:
        if k not in cur: raise SystemExit(f"Missing key in {p}: {'.'.join(keys)}")
        cur = cur[k]
    return cur

def q(s): return "'" + str(s).replace("'", "'\"'\"'") + "'"

units_family   = req(d,'units','family')
dt_fs          = float(req(d,'time','dt_fs'))
t_total_ps     = float(req(d,'time','t_total_ps'))
bx             = req(d,'boundaries','x'); by = req(d,'boundaries','y'); bz = req(d,'boundaries','z')
r0_nm          = float(req(d,'source','r0_nm'))
Se_eV_per_A    = float(req(d,'source','Se_eV_per_A'))
g0             = float(req(d,'ttm','coupling','g0_W_m3K'))
ke_file        = req(d,'ttm','ke_curve_file')
Ce_file        = req(d,'ttm','Ce_curve_file')
grid_file      = req(d,'ttm','grid_file')
src_file       = 'inputs/ttm/source_profile.yaml'  # keep simple; you already ship this
rim            = d.get('forcefield',{}).get('rim_langevin',{})
gamma_ps_inv   = rim.get('gamma_ps_inv'); rim_width_nm = rim.get('width_nm'); rim_target_K = rim.get('target_K')
data_path      = req(d,'structure','data_path')
thermo_every   = int(d.get('thermo',{}).get('thermo_every', 100))
use_zbl        = 1 if d.get('fine',{}) or d.get('fine',None) is None and d.get('forcefield',{}).get('line') is None and d.get('forcefield',{}).get('zbl',None) is None and d.get('forcefield',{}).get('zbl_overlay', True) else 0

dt_ps = dt_fs * 1e-3
nsteps = int(round(t_total_ps / dt_ps))
rng = random.Random(int(time.time()))
seed_vel = rng.randrange(10**6, 10**9); seed_vac = rng.randrange(10**6, 10**9)

print("UNITS="        + q(units_family))
print("DT="           + q(f"{dt_ps:.12g}"))
print("NSTEPS="       + q(str(nsteps)))
print("R0="           + q(f"{r0_nm:.12g}"))
print("SE="           + q(f"{Se_eV_per_A:.12g}"))
print("G0="           + q(f"{g0:.12g}"))
print("BX="           + q(bx)); print("BY=" + q(by)); print("BZ=" + q(bz))
print("DATAFILE="     + q(data_path))
print("KE_FILE="      + q(ke_file))
print("CE_FILE="      + q(Ce_file))
print("GRID_FILE="    + q(grid_file))
print("SRC_FILE="     + q(src_file))
print("THERMO_EVERY=" + q(str(thermo_every)))
print("USE_ZBL="      + q(str(use_zbl)))
print("SEED_VEL="     + q(str(seed_vel)))
print("SEED_VAC="     + q(str(seed_vac)))
print("PHYS_PATH="    + q(str(p)))
if gamma_ps_inv is not None: print("RIM_GAMMA="    + q(f"{float(gamma_ps_inv):.12g}"))
if rim_width_nm  is not None: print("RIM_WIDTH_NM=" + q(f"{float(rim_width_nm):.12g}"))
if rim_target_K  is not None: print("RIM_TK="       + q(f"{float(rim_target_K):.12g}"))
PYCODE
)" || { echo "ERROR: failed to parse $PHYS_YAML" >&2; exit 3; }

# -------- Resolve run dir & snapshot --------
timestamp="$(date +%Y%m%d_%H%M%S)"
: "${RUN_DIR:=runs/phase1_baseline/${timestamp}}"

SNAP="$RUN_DIR/inputs_snapshot"
mkdir -p "$RUN_DIR"/{logs,dumps,post,figures} \
         "$SNAP"/{lammps/includes,structure,ttm}

# Snapshot inputs
cp -f "$PHYS_PATH"                    "$SNAP/physics.yaml"        || true
cp -f config/metrics.yaml             "$SNAP/metrics.yaml"        || true
cp -f config/figure_style.yaml        "$SNAP/figure_style.yaml"   || true
cp -f "$DATAFILE"                     "$SNAP/structure/"          || true
cp -f inputs/lammps/in.main.lmp       "$SNAP/lammps/"
cp -f inputs/lammps/in.equil.lmp      "$SNAP/lammps/"
cp -f inputs/lammps/in.ttm_spike.lmp  "$SNAP/lammps/"
cp -f inputs/lammps/includes/*.in     "$SNAP/lammps/includes/" 2>/dev/null || true
[[ -f inputs/lammps/tables/zbl.table ]] && cp -f inputs/lammps/tables/zbl.table "$SNAP/lammps/"
cp -f "$KE_FILE"                      "$SNAP/ttm/" 2>/dev/null || true
cp -f "$CE_FILE"                      "$SNAP/ttm/" 2>/dev/null || true
cp -f "$GRID_FILE"                    "$SNAP/ttm/" 2>/dev/null || true
cp -f "$SRC_FILE"                     "$SNAP/ttm/" 2>/dev/null || true

# -------- Precompute Te.in (cylindrical Gaussian on TTM grid) --------
TEF="${RUN_DIR}/post/Te.in"
python3 - <<PYGEN || { echo "ERROR: Te.in generation failed" >&2; exit 88; }
import os, math

DATAFILE = r"""${DATAFILE}"""
R0_nm    = float(r"""${R0}""")
TE_BASE  = 300.0
TE_CORE  = 20000.0
grid_A   = 1.0  # target spacing in Angstrom

# parse box bounds from LAMMPS data file
xlo=xhi=ylo=yhi=zlo=zhi=None
with open(DATAFILE,"r") as f:
    for ln in f:
        sp=ln.split()
        if len(sp)==4 and sp[2:]==["xlo","xhi"]: xlo,xhi = float(sp[0]), float(sp[1])
        elif len(sp)==4 and sp[2:]==["ylo","yhi"]: ylo,yhi = float(sp[0]), float(sp[1])
        elif len(sp)==4 and sp[2:]==["zlo","zhi"]: zlo,zhi = float(sp[0]), float(sp[1])
        if None not in (xlo,xhi,ylo,yhi,zlo,zhi): break
if None in (xlo,xhi,ylo,yhi,zlo,zhi):
    raise SystemExit("could not parse box from data file")

Lx, Ly, Lz = xhi-xlo, yhi-ylo, zhi-zlo
NX = int(math.ceil(Lx / grid_A)); NY = int(math.ceil(Ly / grid_A)); NZ = int(math.ceil(Lz / grid_A))
dx, dy = Lx/NX, Ly/NY
sigma = (R0_nm*10.0)/math.sqrt(2.0)  # convert nm->Å then Gaussian sigma

os.makedirs(os.path.dirname(r"""${TEF}"""), exist_ok=True)
with open(r"""${TEF}""","w") as out:
    for iz in range(1, NZ+1):
        for iy in range(1, NY+1):
            y = (iy-0.5)*dy - 0.5*Ly
            for ix in range(1, NX+1):
                x = (ix-0.5)*dx - 0.5*Lx
                r2 = x*x + y*y
                Te = TE_BASE + (TE_CORE-TE_BASE)*math.exp(-r2/(2.0*sigma*sigma))
                out.write(f"{ix} {iy} {iz} {Te:.6f}\n")

with open(r"""${TEF}"""+".meta.txt","w") as m:
    m.write(f"NX NY NZ = {NX} {NY} {NZ}\n")
    m.write(f"Lx Ly Lz (A) = {Lx:.6f} {Ly:.6f} {Lz:.6f}\n")
    m.write(f"grid_A = {grid_A}\n")
    m.write(f"R0_nm = {R0_nm}\n")
PYGEN

# quick sanity
wc -l "$TEF" || true
head -n 2 "$TEF" || true
tail -n 2 "$TEF" || true

# -------- Manifest --------
cat > "$RUN_DIR/manifest.yaml" <<EOF
phase: 1
units: ${UNITS}
dt_ps: ${DT}
nsteps: ${NSTEPS}
thermo_every: ${THERMO_EVERY}
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
  te_infile: ${TEF}
boundaries: { x: ${BX}, y: ${BY}, z: ${BZ} }
seeds:
  velocity: ${SEED_VEL}
  vacancy:  ${SEED_VAC}
snapshot_dir: ${SNAP}
timestamp: ${timestamp}
EOF

# -------- Build LAMMPS command --------
LAMMPS_IN=inputs/lammps/in.main.lmp
LOGFILE="$RUN_DIR/logs/log.lammps"
SCREEN="$RUN_DIR/logs/screen.out"
: > "$SCREEN"; : > "$LOGFILE"; mkdir -p "$(dirname "$SCREEN")" "$(dirname "$LOGFILE")"

LAMMPS_CMD=(
  "$LMP_BIN" -echo both
  -in "$LAMMPS_IN" -log "$LOGFILE"
  -var UNITS "$UNITS" -var DATAFILE "$DATAFILE"
  -var DT "$DT" -var NSTEPS "$NSTEPS"
  -var R0 "$R0" -var SE "$SE" -var G0 "$G0"
  -var KE_FILE "$KE_FILE" -var CE_FILE "$CE_FILE"
  -var GRID_FILE "$GRID_FILE" -var SRC_FILE "$SRC_FILE"
  -var BX "$BX" -var BY "$BY" -var BZ "$BZ"
  -var OUTDIR "$RUN_DIR"
  -var THERMO_EVERY "$THERMO_EVERY"
  -var USE_ZBL "$USE_ZBL"
  -var SEED_VEL "$SEED_VEL"
  -var TEF "$TEF"
)
[[ -n "${RIM_GAMMA-}"    ]] && LAMMPS_CMD+=( -var RIM_GAMMA "$RIM_GAMMA" )
[[ -n "${RIM_WIDTH_NM-}" ]] && LAMMPS_CMD+=( -var RIM_WIDTH_NM "$RIM_WIDTH_NM" )
[[ -n "${RIM_TK-}"       ]] && LAMMPS_CMD+=( -var RIM_TK "$RIM_TK" )

echo "=== Phase 1 baseline ==="
echo "Run dir : $RUN_DIR"
echo "lmp     : $LMP_BIN"
echo "in.lmp  : $LAMMPS_IN"
echo "DT(ps)  : $DT   NSTEPS: $NSTEPS   r0(nm): $R0   Se(eV/Å): $SE   g0: $G0"
echo "thermo_every: $THERMO_EVERY   use ZBL: $USE_ZBL"
echo "Te.in   : $TEF"
echo

if [[ $DRY_RUN -eq 1 ]]; then
  printf '%q ' "${LAMMPS_CMD[@]}"; echo
  echo "(dry-run only)"
  exit 0
fi

# -------- Launch under SLURM (default -n 1 like Phase 0) --------
LAUNCH=()
if command -v srun >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
  LAUNCH=(srun -u -n "${SLURM_NTASKS:-1}")
fi

STD_BUF=()
if command -v stdbuf >/dev/null 2>&1; then
  STD_BUF=(stdbuf -oL -eL)
fi

"${LAUNCH[@]}" "${STD_BUF[@]}" "${LAMMPS_CMD[@]}" >"$SCREEN" 2>&1 || {
  rc=$?
  echo "LAMMPS exited with RC=$rc" >&2
  echo "---- screen.out (tail) ----" >&2; tail -n 120 "$SCREEN" >&2 || true
  echo "---- log.lammps (tail) ----" >&2;  tail -n 120 "$LOGFILE" >&2 || true
  exit "$rc"
}

echo "Done. Logs at: $RUN_DIR/logs"
echo "Snapshot at:   $SNAP"
