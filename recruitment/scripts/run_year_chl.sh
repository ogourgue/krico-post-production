#!/bin/bash
#SBATCH --job-name=chlorophyll
#SBATCH --qos=nf
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --array=1-32
#SBATCH --output=logs/chlorophyll_%A_%a.out
#SBATCH --error=logs/chlorophyll_%A_%a.err

# ============================================================================
# KRICO chlorophyll indices, pass 2 of dataset v3.0.0 — SLURM array
# ============================================================================
#
# One array task per spawning year. Spawning years 1994-2025 (32 years) map
# to array indices 1-32, and each year spans 4 simulation folders, as in
# run_year.sh.
#
# Within each folder, each daily release file YYYY_MM_DD.nc goes in series
# through process_cohort_chl.py. The output folder is named on the command
# line, relative to recruitment/, and must exist; its name must end in
# _chl (or _chl_test), so that pass 2 can never write into the folder of a
# classification.
#
# Besides the trajectories, a cohort needs the satellite climatology
# (download_globcolour.sh) and the monthly files of the model
# (download_pisces.sh). Neither the pass 1 files nor the release-cell files
# are read.
#
# Existing output files are skipped, so a failed or interrupted task is
# simply resubmitted. A file is written under a temporary name and renamed
# when complete: a file with its final name is a complete file.
#
# The files record the commit that produced them: a cohort stops if tracked
# files are modified or if the commit is not pushed.
#
# Before submitting, check the inputs of every cohort (check_inputs_chl.sh).
#
# Submit from recruitment/scripts:
#   mkdir -p logs ../data_v3.0_chl
#   sbatch --array=1 run_year_chl.sh data_v3.0_chl     # 1994, as a test
#   sbatch --array=2-32 run_year_chl.sh data_v3.0_chl  # the rest
#
# Requires KRICO_RUNS, KRICO_GLOBCOLOUR and KRICO_PISCES.
# ============================================================================

set -euo pipefail

# ----------------------------------------------------------------------------
# Environment
# ----------------------------------------------------------------------------
module load python3

for name in KRICO_RUNS KRICO_GLOBCOLOUR KRICO_PISCES; do
    if [[ -z "${!name:-}" ]]; then
        echo "ERROR: ${name} environment variable not set." >&2
        exit 1
    fi
done

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
# Use SLURM_SUBMIT_DIR when running under SLURM (BASH_SOURCE does not resolve
# correctly when sbatch copies the script to the SLURM spool directory).
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
RECRUITMENT_DIR="$(dirname "$SCRIPT_DIR")"
SCRIPT="${SCRIPT_DIR}/process_cohort_chl.py"

FOLDER="${1:?give the output folder under recruitment/, e.g. data_v3.0_chl}"
DATA_DIR="${RECRUITMENT_DIR}/${FOLDER}"
if [[ ! -d "${DATA_DIR}" ]]; then
    echo "ERROR: output folder ${DATA_DIR} does not exist; create it first." >&2
    exit 1
fi

# Pass 2 can never write into the folder of a classification.
if [[ "${FOLDER}" != *"_chl" && "${FOLDER}" != *"_chl_test" ]]; then
    echo "ERROR: the output folder of pass 2 must end in _chl or _chl_test." >&2
    exit 1
fi

# ----------------------------------------------------------------------------
# Resolve spawning year and the four simulation folders that feed it
# ----------------------------------------------------------------------------
YEAR_INDEX=${SLURM_ARRAY_TASK_ID:?run through sbatch --array}   # 1..32
SPAWNING_YEAR=$((1993 + YEAR_INDEX))                            # 1994..2025

FIRST_FOLDER=$((4 * (SPAWNING_YEAR - 1994) + 1))
FOLDERS=()
for offset in 0 1 2 3; do
    n=$((FIRST_FOLDER + offset))
    FOLDERS+=("$(printf 'KRICO_%04d' "$n")")
done

echo "================================================================"
echo "Chlorophyll indices (pass 2) for spawning year ${SPAWNING_YEAR}"
echo "Array task: ${YEAR_INDEX}/32"
echo "Folders: ${FOLDERS[*]}"
echo "Trajectory input: ${KRICO_RUNS}"
echo "Satellite climatology: ${KRICO_GLOBCOLOUR}"
echo "Model hindcast: ${KRICO_PISCES}"
echo "Output: ${DATA_DIR}"
echo "Script: ${SCRIPT}"
echo "================================================================"

# ----------------------------------------------------------------------------
# Chlorophyll indices of all cohorts of this spawning year
# ----------------------------------------------------------------------------
n_processed=0
n_skipped=0
n_failed=0

for folder in "${FOLDERS[@]}"; do
    folder_path="${KRICO_RUNS}/${folder}"

    if [[ ! -d "${folder_path}" ]]; then
        echo "WARNING: folder not found, skipping: ${folder_path}" >&2
        continue
    fi

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
        if python3 -u "${SCRIPT}" "${cohort_file}" "${output_file}"; then
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
