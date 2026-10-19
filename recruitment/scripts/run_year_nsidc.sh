#!/bin/bash
#SBATCH --job-name=recruitment_nsidc
#SBATCH --qos=nf
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --array=1-32
#SBATCH --output=logs/recruitment_nsidc_%A_%a.out
#SBATCH --error=logs/recruitment_nsidc_%A_%a.err

# ============================================================================
# KRICO recruitment post-processing with NSIDC sea ice — SLURM array
# ============================================================================
#
# Same classification as run_year.sh, but every sea-ice concentration it
# reads (M1 at spawning, and the along-trajectory series behind M4, the
# advance detection and censoring) comes from the NSIDC G02202 v6 record
# instead of GLORYS12, sampled at the trajectory positions and days
# (krico_recruitment.nsidc_ice). The trajectories, and so the circulation,
# are unchanged. This is the ice-source test for the ERA-Interim/ERA5
# forcing change in GLORYS12 at the end of 2018. The ice source is fixed in
# this script rather than passed in, so that the file itself records what
# was run. Output goes to recruitment/data_nsidc/, never to data/.
#
# Coast and land cells read NaN (they are fill in the CDR); the nearest-ocean
# redirection of process_cohort.py is deliberately not used, after
# validate_nsidc.py showed coast cells hold 0.4% of positions and move M1 by
# at most 0.4 points (validate_nsidc.csv).
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
# series. Outputs go to recruitment/data_nsidc/ as YYYY_MM_DD.nc (release date in
# filename uniquely identifies the cohort across folders). Job logs go to
# recruitment/scripts/logs/, since the SBATCH output paths above are relative
# to the submission directory.
#
# Existing output files are skipped (idempotent re-runs). Note that this makes
# the script a no-op against a populated output directory: when the
# classification itself changes, move the previous recruitment/data_nsidc/
# aside first so that a fresh one is created.
#
# M1 is evaluated at spawning, which precedes release by the descent-ascent
# interval (krico_recruitment.sea_ice.SPAWNING_OFFSET_DAYS, echoed by
# process_cohort.py for each cohort), from the NSIDC monthly stacks
# nsidc_ice_YYYY_MM.nc written by preprocess_nsidc.py.
#
# Submit from recruitment/scripts:
#   export KRICO_RUNS=/scratch/cvan/KRICO/Runs
#   export KRICO_NSIDC=/scratch/cvan/KRICO/Pre/NSIDC
#   sbatch run_year_nsidc.sh
# Or for a subset (e.g. just first year):
#   sbatch --array=1 run_year_nsidc.sh
# ============================================================================

set -euo pipefail

# ----------------------------------------------------------------------------
# Environment
# ----------------------------------------------------------------------------
module load python3

# Validate KRICO_RUNS environment variable
if [[ -z "${KRICO_RUNS:-}" ]]; then
    echo "ERROR: KRICO_RUNS environment variable not set." >&2
    echo "Set it to the directory containing raw trajectory simulations, e.g.:" >&2
    echo "  export KRICO_RUNS=/scratch/cvan/KRICO/Runs" >&2
    exit 1
fi

# Validate KRICO_NSIDC environment variable: the monthly NSIDC stacks and
# the ancillary/ grid file, both needed by the sampler.
if [[ -z "${KRICO_NSIDC:-}" ]]; then
    echo "ERROR: KRICO_NSIDC environment variable not set." >&2
    echo "Set it to the NSIDC preprocessing directory containing" >&2
    echo "nsidc_ice_YYYY_MM.nc and ancillary/, e.g.:" >&2
    echo "  export KRICO_NSIDC=/scratch/cvan/KRICO/Pre/NSIDC" >&2
    exit 1
fi

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
# Fixed directory structure: script is in recruitment/scripts/
# Use SLURM_SUBMIT_DIR when running under SLURM (BASH_SOURCE does not resolve
# correctly when sbatch copies the script to the SLURM spool directory)
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
RECRUITMENT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="${RECRUITMENT_DIR}/data_nsidc"
SCRIPT="${SCRIPT_DIR}/process_cohort.py"

mkdir -p "${DATA_DIR}"

# SLURM resolves the SBATCH output paths relative to the submission directory,
# so the log directory must be created there rather than in RECRUITMENT_DIR.
mkdir -p "${SCRIPT_DIR}/logs"

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
echo "Recruitment post-processing (NSIDC sea ice) for spawning year ${SPAWNING_YEAR}"
echo "Array task: ${YEAR_INDEX}/32"
echo "Folders: ${FOLDERS[*]}"
echo "Trajectory input: ${KRICO_RUNS}"
echo "Sea-ice input: ${KRICO_NSIDC} (NSIDC G02202 v6)"
echo "Recruitment output: ${DATA_DIR}"
echo "Script: ${SCRIPT}"
echo "================================================================"

# ----------------------------------------------------------------------------
# Process all cohorts for this spawning year
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
        if python3 "${SCRIPT}" "${cohort_file}" "${output_file}" \
                   --ice-source nsidc --nsidc "${KRICO_NSIDC}"; then
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
