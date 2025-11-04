#!/bin/bash
#SBATCH -J ceo2_p1_baseline
#SBATCH -p eertekin
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -D /u/gclopton/scratch/lammps
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err


set -Eeuo pipefail
IFS=$'\n\t'
trap 'echo "[ERR] line $LINENO: $BASH_COMMAND" >&2' ERR

# ---- Toolchain: SAME as Phase 0 (lmp was built with GCC/OpenMPI/FFTW) ----
module purge
module load gcc/13.3.0
module load openmpi/5.0.1
module load fftw/3.3.10

# ---- Python env (PyYAML & plotting) ------
source /sw/apps/anaconda3/2024.10/etc/profile.d/conda.sh
conda activate ceo2-shi

# ------------ Paths & binaries ----
REPO_DIR="$PWD"
RUNNER="${REPO_DIR}/scripts/run_nominal.sh"
LAMMPS_ROOT="${REPO_DIR}/opt/lammps/stable_29Aug2024_update1"
LMP_BIN="${LAMMPS_ROOT}/bin/lmp"
PHYS="${REPO_DIR}/config/physics.yaml"

export LD_LIBRARY_PATH="${LAMMPS_ROOT}/lib64:${LAMMPS_ROOT}/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

# ---------- Dated run folder under runs/phase1_baseline ----
STAMP="$(date +%b_%d_%H%M%S)"
RUN_TAG="${STAMP}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
RUN_DIR="${REPO_DIR}/runs/phase1_baseline/${RUN_TAG}"
mkdir -p "${RUN_DIR}/logs" "${RUN_DIR}/dumps" "${RUN_DIR}/post" "${RUN_DIR}/figures" "${RUN_DIR}/inputs_snapshot"

# Log this job’s stdout/err into the run folder (like Phase 0)
exec > >(tee -a "${RUN_DIR}/logs/job.out") 2> >(tee -a "${RUN_DIR}/logs/job.err" >&2)
echo "Run tag : ${RUN_TAG}"
echo "Node(s): ${SLURM_NODELIST}"
echo "lmp    : ${LMP_BIN}"
ls -l "${RUNNER}" || { echo "ERROR: runner not found or not executable"; exit 127; }

# ============ Snapshot inputs for provenance ===========
rsync -a --delete "${REPO_DIR}/config/"            "${RUN_DIR}/inputs_snapshot/config/"
rsync -a --delete "${REPO_DIR}/inputs/lammps/"     "${RUN_DIR}/inputs_snapshot/lammps/"
rsync -a --delete "${REPO_DIR}/inputs/ttm/"        "${RUN_DIR}/inputs_snapshot/ttm/"
rsync -a --delete "${REPO_DIR}/inputs/structure/"  "${RUN_DIR}/inputs_snapshot/structure/"

# ==================== Launch Phase 1 via the SAME runner as Phase 0 ===============
bash -x "${RUNNER}" --lmp "${LMP_BIN}" --config "${PHYS}" --run-dir "${RUN_DIR}"
RC=$?

# ---- Post-processing (Phase-1 script if present, else Phase-0) ---- #TODO: Remove Phase 0
if [[ -x "${REPO_DIR}/scripts/postprocess_phase1.sh" ]]; then
  RUN="${RUN_DIR}" "${REPO_DIR}/scripts/postprocess_phase1.sh" || true
elif [[ -x "${REPO_DIR}/scripts/postprocess_phase0.sh" ]]; then
  RUN="${RUN_DIR}" "${REPO_DIR}/scripts/postprocess_phase0.sh" || true
fi

# Prefer the energy TSV if present
THERMO_TSV="${RUN_DIR}/post/thermo_energy.tsv"
if [[ -f "${THERMO_TSV}" ]]; then
  python "${REPO_DIR}/scripts/energy_conditional.py" >/dev/null 2>&1 || true  # no-op if absent
  python "${REPO_DIR}/scripts/energy_balance.py" \
         --run "${RUN_DIR}" --thermo "${THERMO_TSV}" \
         --physics "${PHYS}" --tol 0.02 || true
else
  python "${REPO_DIR}/scripts/energy_balance.py" \
         --run "${RUN_DIR}" \
         --physics "${PHYS}" --tol 0.02 || true
fi

# ---- LATEST symlink and exit with runner’s status -----------------
mkdir -p "${REPO_DIR}/runs/phase1_baseline"
ln -sfn "${RUN_DIR}" "${REPO_DIR}/runs/phase1_baseline/LATEST"
exit "${RC}"
