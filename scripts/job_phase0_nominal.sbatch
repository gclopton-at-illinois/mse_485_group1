#!/bin/bash
#SBATCH -J ceo2_phase0
#SBATCH -p mrsec
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
# Always start in the repo root so relative paths resolve consistently
#SBATCH -D /u/gclopton/scratch/lammps
# Write SLURM’s logs to flat files (created by SLURM at launch)
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
# #SBATCH -A <your_allocation>   # uncomment & set if SLURM requires an account

set -Eeuo pipefail
IFS=$'\n\t'
trap 'echo "[ERR] line $LINENO: $BASH_COMMAND" >&2' ERR

# ---------- Resolve repo root & runner (honor #SBATCH -D via $PWD) ----------
REPO_DIR="$PWD"
RUNNER="${REPO_DIR}/scripts/run_nominal.sh"

echo "SUBMIT_DIR: ${SLURM_SUBMIT_DIR:-<unset>}"
echo "PWD:        ${PWD}"
echo "REPO_DIR:   ${REPO_DIR}"
echo "Runner:     ${RUNNER}"
ls -l "${RUNNER}" || { echo "ERROR: runner not found or not executable"; exit 127; }

# ---------- Toolchain (match your build) ----------
module purge
module load gcc/13.3.0
module load openmpi/5.0.1
module load fftw/3.3.10

# ---------- Python / Conda env (needed for PyYAML in the runner) ----------
source /sw/apps/anaconda3/2024.10/etc/profile.d/conda.sh
conda activate ceo2-shi

# Prove PyYAML is present (prints into logs)
python - <<'PY'
import sys, yaml
print(f"Python: {sys.version.split()[0]}")
print(f"PyYAML: {yaml.__version__}")
PY

# ---------- LAMMPS (shared) ----------
LAMMPS_ROOT=/u/gclopton/scratch/lammps/opt/lammps/stable_29Aug2024_update1
LMP_BIN="${LAMMPS_ROOT}/bin/lmp"
export LD_LIBRARY_PATH="${LAMMPS_ROOT}/lib64:${LAMMPS_ROOT}/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

# ---------- Dated run folder (always) ----------
STAMP="$(date +%b_%d_%H%M%S)"                         # e.g., Oct_23_143855
RUN_TAG="${STAMP}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}"  # e.g., Oct_23_143855_ceo2_phase0_5492143
RUN_DIR="${REPO_DIR}/runs/phase0_nominal/${RUN_TAG}"
mkdir -p "${RUN_DIR}/logs" "${RUN_DIR}/dumps" "${RUN_DIR}/post"

# Redirect ALL script output into the dated folder
exec > >(tee -a "${RUN_DIR}/logs/job.out") 2> >(tee -a "${RUN_DIR}/logs/job.err" >&2)
echo "Run tag: ${RUN_TAG}"
echo "Nodes:   ${SLURM_NODELIST}"
echo "LAMMPS:  ${LMP_BIN}"

# ---------- Launch Phase 0 via runner ----------
# No outer srun here; the runner will do the single-task srun for lmp.
set +e
bash -x "${RUNNER}" --lmp "${LMP_BIN}" --run-dir "${RUN_DIR}"
RC=$?
set -e

# ---------- Mirror SLURM’s flat logs into the dated folder ----------
# (These files are created by SLURM at job start thanks to -o/-e above.)
SLURM_OUT="${REPO_DIR}/slurm-${SLURM_JOB_ID}.out"
SLURM_ERR="${REPO_DIR}/slurm-${SLURM_JOB_ID}.err"
[[ -f "${SLURM_OUT}" ]] && cp -f "${SLURM_OUT}" "${RUN_DIR}/logs/slurm.out" || true
[[ -f "${SLURM_ERR}" ]] && cp -f "${SLURM_ERR}" "${RUN_DIR}/logs/slurm.err" || true

# ---------- Convenience symlink ----------
ln -sfn "${RUN_DIR}" "${REPO_DIR}/runs/phase0_nominal/LATEST"

echo "Done. See ${RUN_DIR}/logs and ${RUN_DIR}/dumps."
exit "${RC}"
