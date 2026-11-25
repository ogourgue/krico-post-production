#!/bin/bash
#SBATCH --job-name=recruitment
#SBATCH --qos=nf
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --array=1-32
#SBATCH --output=logs/recruitment_%A_%a.out
#SBATCH --error=logs/recruitment_%A_%a.err

# ============================================================================
# KRICO recruitment classification, pass 1 of dataset v3.0.0 — SLURM array
# ============================================================================
#
# One array task per spawning year. Spawning years 1994-2025 (32 years) map
# to array indices 1-32.
#
# Each spawning year spans 4 simulation folders (releases of 15 November to
# 14 March):
#   KRICO_{4(Y-1994)+1} to KRICO_{4(Y-1994)+4}
#
# Within each folder, each daily release file YYYY_MM_DD.nc is classified in
# series by process_cohort.py. The output folder is named on the command
# line, relative to recruitment/, and must exist: there is no default, so
# that a classification can never land in the folder of another one.
#
# Besides the trajectories, a cohort needs the release-cell file of its
# spawning year (preprocess_glorys12_release_cells.sh), from which M1, M3
# and the temperature at spawning are read. No GLORYS12 forcing file is
# read.
#
# Existing output files are skipped, so a failed or interrupted task is
# simply resubmitted. A file is written under a temporary name and renamed
# when complete: a file with its final name is a complete file.
#
# The files record the commit that produced them: a cohort stops if tracked
# files are modified or if the commit is not pushed.
#
# Before submitting, check the inputs of every cohort (check_inputs.sh).
#
# Submit from recruitment/scripts:
#   mkdir -p logs ../data_v3.0_glorys12
#   sbatch --array=1 run_year.sh data_v3.0_glorys12     # 1994, as a test
#   sbatch --array=2-32 run_year.sh data_v3.0_glorys12  # the rest
#
# A second word names another ice source (the ice-source test: see
# process_cohort.py); each source has its own folder:
#   mkdir -p ../data_v3.0_nsidc ../data_v3.0_glorys12_25km
#   sbatch run_year.sh data_v3.0_nsidc nsidc
#   sbatch run_year.sh data_v3.0_glorys12_25km glorys12_25km
#
# Requires KRICO_RUNS and KRICO_GLORYS12_RELEASE_CELLS, and KRICO_NSIDC or
# KRICO_GLORYS12_25KM for the other ice sources.
# ============================================================================

set -euo pipefail

# ----------------------------------------------------------------------------
# Environment
# ----------------------------------------------------------------------------
module load python3

for name in KRICO_RUNS KRICO_GLORYS12_RELEASE_CELLS; do
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
SCRIPT="${SCRIPT_DIR}/process_cohort.py"

FOLDER="${1:?give the output folder under recruitment/, e.g. data_v3.0_glorys12}"
DATA_DIR="${RECRUITMENT_DIR}/${FOLDER}"
if [[ ! -d "${DATA_DIR}" ]]; then
    echo "ERROR: output folder ${DATA_DIR} does not exist; create it first." >&2
    exit 1
fi

# Ice source, and the folder that must go with it: a classification with
# one sea ice can never be written into the folder of another.
ICE_SOURCE="${2:-glorys12}"
case "${ICE_SOURCE}" in
    glorys12)       ICE_ENV="" ;;
    nsidc)          ICE_ENV="KRICO_NSIDC" ;;
    glorys12_25km)  ICE_ENV="KRICO_GLORYS12_25KM" ;;
    *)
        echo "ERROR: unknown ice source '${ICE_SOURCE}' (glorys12, nsidc or glorys12_25km)" >&2
        exit 1
        ;;
esac
if [[ "${FOLDER}" != *"_${ICE_SOURCE}" ]]; then
    echo "ERROR: the output folder must end in _${ICE_SOURCE} for ice source ${ICE_SOURCE}." >&2
    exit 1
fi
if [[ -n "${ICE_ENV}" && -z "${!ICE_ENV:-}" ]]; then
    echo "ERROR: ${ICE_ENV} environment variable not set." >&2
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
echo "Recruitment classification (pass 1) for spawning year ${SPAWNING_YEAR}"
echo "Array task: ${YEAR_INDEX}/32"
echo "Folders: ${FOLDERS[*]}"
echo "Trajectory input: ${KRICO_RUNS}"
echo "Release-cell input: ${KRICO_GLORYS12_RELEASE_CELLS}"
echo "Sea ice: ${ICE_SOURCE}${ICE_ENV:+ (${!ICE_ENV})}"
echo "Output: ${DATA_DIR}"
echo "Script: ${SCRIPT}"
echo "================================================================"

# ----------------------------------------------------------------------------
# Classify all cohorts of this spawning year
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
        if python3 -u "${SCRIPT}" "${cohort_file}" "${output_file}" \
                   --ice-source "${ICE_SOURCE}"; then
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
