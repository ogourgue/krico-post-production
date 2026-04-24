#!/bin/bash
#SBATCH --job-name=recruitment
#SBATCH --qos=nf
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --array=1-32
#SBATCH --output=logs/recruitment_%A_%a.out
#SBATCH --error=logs/recruitment_%A_%a.err

# ============================================================================
# KRICO recruitment post-processing — SLURM array
# ============================================================================
#
# One array task per spawning year. Spawning years 1994-2025 (32 years) map
# to array indices 1-32.
#
# Each spawning year spans 4 simulation folders (Nov, Dec, Jan, Feb releases):
#   KRICO_{4(Y-1994)+1}  Nov releases (year Y-1)
#   KRICO_{4(Y-1994)+2}  Dec releases (year Y-1)
#   KRICO_{4(Y-1994)+3}  Jan releases (year Y)
#   KRICO_{4(Y-1994)+4}  Feb releases (year Y)
#
# Within each folder, each daily release file YYYY_MM_DD.nc is processed in
# series. Outputs go to recruitment/data/ as YYYY_MM_DD.nc (release date in
# filename uniquely identifies the cohort across folders).
#
# Existing output files are skipped (idempotent re-runs).
#
# Submit with:
#   sbatch run_year.sh
# Or for a subset (e.g. just first year):
#   sbatch --array=1 run_year.sh
# ============================================================================

set -euo pipefail

# ----------------------------------------------------------------------------
# Environment
# ----------------------------------------------------------------------------
module load python3
export KRICO_ROOT=/scratch/cvan/KRICO

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
if [[ -z "${KRICO_ROOT:-}" ]]; then
    echo "ERROR: KRICO_ROOT environment variable not set." >&2
    exit 1
fi

RUNS_DIR="${KRICO_ROOT}/Runs"
RECRUITMENT_DIR="${KRICO_ROOT}/Post/Production/recruitment"
DATA_DIR="${RECRUITMENT_DIR}/data"
SCRIPT="${RECRUITMENT_DIR}/scripts/process_cohort.py"

mkdir -p "${DATA_DIR}"
mkdir -p "${RECRUITMENT_DIR}/logs"

# ----------------------------------------------------------------------------
# Resolve spawning year and the four simulation folders that feed it
# ----------------------------------------------------------------------------
YEAR_INDEX=${SLURM_ARRAY_TASK_ID:-1}            # 1..32
SPAWNING_YEAR=$((1993 + YEAR_INDEX))            # 1994..2025

# KRICO_NNNN folders: 4 per spawning year, indexed sequentially.
FIRST_FOLDER=$((4 * (SPAWNING_YEAR - 1994) + 1))
FOLDERS=()
for offset in 0 1 2 3; do
    n=$((FIRST_FOLDER + offset))
    FOLDERS+=("$(printf 'KRICO_%04d' "$n")")
done

echo "================================================================"
echo "Recruitment post-processing for spawning year ${SPAWNING_YEAR}"
echo "Array task: ${YEAR_INDEX}/32"
echo "Folders: ${FOLDERS[*]}"
echo "Data output: ${DATA_DIR}"
echo "Script: ${SCRIPT}"
echo "================================================================"

# ----------------------------------------------------------------------------
# Process all cohorts for this spawning year
# ----------------------------------------------------------------------------
n_processed=0
n_skipped=0
n_failed=0

for folder in "${FOLDERS[@]}"; do
    folder_path="${RUNS_DIR}/${folder}"

    if [[ ! -d "${folder_path}" ]]; then
        echo "WARNING: folder not found, skipping: ${folder_path}" >&2
        continue
    fi

    # Glob all daily release files in this folder.
    shopt -s nullglob
    cohort_files=("${folder_path}"/[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9].nc)
    shopt -u nullglob

    if [[ ${#cohort_files[@]} -eq 0 ]]; then
        echo "WARNING: no cohort files in ${folder_path}" >&2
        continue
    fi

    echo
    echo "--- Folder: ${folder} (${#cohort_files[@]} cohorts) ---"

    for cohort_file in "${cohort_files[@]}"; do
        basename="$(basename "${cohort_file}")"          # e.g. 1993_11_15.nc
        output_file="${DATA_DIR}/${basename}"

        if [[ -f "${output_file}" ]]; then
            echo "  [skip] ${basename} (output exists)"
            n_skipped=$((n_skipped + 1))
            continue
        fi

        echo "  [run]  ${basename}"
        if python3 "${SCRIPT}" "${cohort_file}" "${output_file}"; then
            n_processed=$((n_processed + 1))
        else
            echo "  [FAIL] ${basename}" >&2
            n_failed=$((n_failed + 1))
            # Continue with the next cohort rather than aborting the whole year.
        fi
    done
done

echo
echo "================================================================"
echo "Spawning year ${SPAWNING_YEAR} complete"
echo "  processed: ${n_processed}"
echo "  skipped:   ${n_skipped} (output already existed)"
echo "  failed:    ${n_failed}"
echo "================================================================"

if [[ ${n_failed} -gt 0 ]]; then
    exit 1
fi
