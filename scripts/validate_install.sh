#!/usr/bin/env bash
# LAMMPS feature & smoke-test validator (robust to stderr/help format)
# Usage:
#   scripts/validate_install.sh [/full/path/to/lmp]
#
# This script now attempts to load the same module stack used by jobs
# (gcc/13.3.0, openmpi/5.0.1, fftw/3.3.10) automatically if Lmod is
# available. Set VALIDATOR_LOAD_MODULES=0 to disable, or override the
# module list with VALIDATOR_MODULES.
#
# Exits nonzero on failure; prints brief diagnostics.

set -Eeuo pipefail

# --- optional: load job-like modules via Lmod ---
VALIDATOR_LOAD_MODULES=${VALIDATOR_LOAD_MODULES:-1}
VALIDATOR_MODULES=${VALIDATOR_MODULES:-"gcc/13.3.0 openmpi/5.0.1 fftw/3.3.10"}

if [[ "$VALIDATOR_LOAD_MODULES" == "1" ]]; then
  # Ensure 'module' is defined (Lmod)
  if ! type module >/dev/null 2>&1; then
    # Try standard Lmod init locations
    for init in \
      /usr/share/lmod/lmod/init/bash \
      /etc/profile.d/lmod.sh \
      /etc/profile.d/modules.sh; do
      [[ -r "$init" ]] && source "$init" && break || true
    done
  fi
  if type module >/dev/null 2>&1; then
    # Silence Lmod debug if enabled in this shell
    unset LMOD_SH_DBG_ON || true
    module purge || true
    # Load declared module stack; ignore missing ones to keep validator portable
    for m in $VALIDATOR_MODULES; do
      module load "$m" || true
    done
  fi
fi

# --- locate lmp ---
if [[ -n "${1-}" ]]; then
  LMP="$1"
else
  LMP="$(command -v lmp || command -v lmp_mpi || command -v lmp_kokkos || command -v lmp_serial || true)"
fi
[[ -x "${LMP:-}" ]] || { echo "ERROR: supply /full/path/to lmp or ensure it is in \$PATH"; exit 1; }

# --- minimal runtime env (libs + quiet MPI on login nodes) ---
bdir="$(dirname "$LMP")"
pfx="$(cd "$bdir/.." && pwd)"   # install prefix that contains bin/, lib/, lib64/
export LD_LIBRARY_PATH="$pfx/lib:$pfx/lib64:${CONDA_PREFIX:+$CONDA_PREFIX/lib}:${LD_LIBRARY_PATH:-}"
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp

echo "Validating: $LMP"
"$LMP" -h | head -n 5 || true

# --- 1) Header capability checks (robust) ---
HELP_FILE="$(mktemp)"; trap 'rm -f "$HELP_FILE"' EXIT
"$LMP" -h 2>&1 >"$HELP_FILE" || true
# Detect dynamic linker/library errors upfront to avoid misleading greps
if grep -qi 'error while loading shared libraries' "$HELP_FILE"; then
  echo "ERROR: $LMP failed to start due to missing shared libraries:" >&2
  cat "$HELP_FILE" >&2
  exit 11
fi

echo -n "[hdr] EXTRA-FIX present? "  ; grep -q "EXTRA-FIX"  "$HELP_FILE" && echo "yes" || { echo "NO"; exit 2; }
echo -n "[hdr] KSPACE present?    "  ; grep -q "KSPACE"     "$HELP_FILE" && echo "yes" || { echo "NO"; exit 3; }
echo -n "[hdr] EXTRA-DUMP present?"; grep -q "EXTRA-DUMP" "$HELP_FILE" && echo "yes" || { echo "NO"; exit 4; }

echo -n "[hdr] TTM fix styles?    "
grep -Eq '^[[:space:]]*ttm($|/grid|/mod)[[:space:]]' "$HELP_FILE" \
  && echo "yes" || { echo "NO"; miss=1; }

# Explicitly assert presence of fix python/invoke (often needed for scripted sources)
echo -n "[hdr] fix python/invoke?  "
grep -Eq '^[[:space:]]*python/invoke([[:space:]]|$)' "$HELP_FILE" \
  && echo "yes" || { echo "NO"; exit 5; }

echo -n "[hdr] thermo_style yaml? "
grep -Eq '^\* Thermo styles|^[[:space:]]*yaml([[:space:]]|$)' "$HELP_FILE" \
  && echo "yes" || { echo "NO"; miss=1; }

echo -n "[hdr] dump ... yaml?     "
grep -Eq '^\* Dump styles|^[[:space:]]*yaml([[:space:]]|$)' "$HELP_FILE" \
  && echo "yes" || { echo "NO"; miss=1; }


# --- 2) Runtime parse/init smoke: PPPM + buck/coul/long + fix ttm (run 0) ---
TMPDIR="$(mktemp -d)"; trap 'rm -rf "$TMPDIR" "$HELP_FILE"' EXIT
TMPIN="$TMPDIR/pppm_ttm.in"; TMPLOG="$TMPDIR/pppm_ttm.log"

cat >"$TMPIN" <<'EOF'
clear
units       metal
atom_style  charge
boundary    p p p
region      box block 0 50 0 50 0 50
create_box  1 box
create_atoms 1 single 10 10 10
create_atoms 1 single 40 40 40
mass 1 140.116
set atom 1 charge  1.0
set atom 2 charge -1.0
neighbor 2.0 bin
kspace_style pppm 1.0e-4
pair_style   buck/coul/long 6.0
pair_coeff   1 1 1000.0 0.5 100.0
group all type 1
# Minimal TTM just to exercise parser/init; numbers are placeholders
fix t all ttm 12345 1.0 1.0 10.0 0.1 0.0 2.0 1 1 1 set 300.0
run 0
EOF

echo "[run] PPPM + buck/coul/long + fix ttm (run 0)…"
if ! "$LMP" -echo none -log "$TMPLOG" -screen none -in "$TMPIN" >/dev/null 2>&1; then
  echo "ERROR: parse/init failed. Last 80 lines:"; tail -n 80 "$TMPLOG" || true; exit 8
fi
echo "      ok"

# --- 3) Runtime YAML smoke: write dump.yaml with embedded thermo ---
YIN="$TMPDIR/yaml_smoke.in"
cat >"$YIN" <<'EOF'
units       lj
atom_style  atomic
lattice     fcc 0.8442
region      box block 0 4 0 4 0 4
create_box  1 box
create_atoms 1 box
mass        1 1.0
pair_style  lj/cut 2.5
pair_coeff  * * 1.0 1.0 2.5
velocity    all create 1.0 4928459
fix         1 all nve
timestep    0.005
thermo      10
thermo_style yaml
dump        d all yaml 20 dump.yaml id type x y z
dump_modify d thermo yes time yes units yes
run         40
EOF

echo "[run] YAML thermo+dump smoke…"
# run LAMMPS *inside* TMPDIR so dump.yaml goes here
(
  cd "$TMPDIR"
  "$LMP" -echo screen -log none -in "$(basename "$YIN")" > screen.out 2>&1
)

if ! grep -q '^thermo:' "$TMPDIR/dump.yaml"; then
  echo "ERROR: dump.yaml written but no embedded 'thermo:' block (missing 'dump_modify d thermo yes'?)"
  echo "Tail of run output:"
  tail -n 80 "$TMPDIR/screen.out" || true
  exit 10
fi
echo "      ok"

echo
echo "All checks passed ✅"
