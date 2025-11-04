#!/bin/bash
#SBATCH -p eertekin                  
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 01:00:00
#SBATCH -J build_lammps_yaml
#SBATCH -o logs/build-%j.out
#SBATCH -e logs/build-%j.err

set -Eeuo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Make sure these exist so Slurm & the build don’t choke
mkdir -p logs .thirdparty opt "$HOME/tmp" "$HOME/.conda/pkgs"

# Campus Cluster often locks /tmp; use $HOME/tmp for all temp work
export TMPDIR="$HOME/tmp"
export TEMP="$HOME/tmp"
export TMP="$HOME/tmp"
export CONDA_PKGS_DIRS="$HOME/.conda/pkgs"

# --- Toolchain modules (CPU build) ---
module purge || true
module load gcc/13.3.0 || module load gcc || true
module load openmpi/5.0.1-gcc-13.3.0 || module load openmpi/5.0.1 || module load openmpi || true
module load cmake || true
module load fftw3 || module load fftw || true

# --- Optional: load/activate Conda env so 'make install-python' can drop the wheel into it ---
module load miniconda3/24.9.2 || true
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate ceo2-shi || echo "[info] continuing without Python env"
fi

# --- Build knobs (override at submit time with: sbatch --export=ALL,TAG=...,WANT_GPU=1 …) ---
export TAG="${TAG:-stable_29Aug2024_update1}"       # a modern stable tag with YAML thermo
export PREFIX="${PREFIX:-$PWD/opt/lammps/${TAG}}"
export SRCROOT="${SRCROOT:-$PWD/.thirdparty}"
export JOBS="${JOBS:-4}"
export WANT_GPU="${WANT_GPU:-0}"                     # 0=CPU (default), 1=KOKKOS/CUDA (only if nodes support it)

bash scripts/build_lammps_yaml.sh
