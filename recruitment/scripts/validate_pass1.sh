#!/bin/bash
#SBATCH --job-name=p1_val
#SBATCH --qos=nf
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/p1_val_%j.out
#SBATCH --error=logs/p1_val_%j.err

# ============================================================================
# Tests of pass 1 on real cohorts — SLURM driver
# ============================================================================
#
# The test is named by a plain word, because ECMWF's sbatch refuses script
# arguments that start with a dash. Submit from recruitment/scripts:
#
#   # regression against v2.0.0 (recruitment/data/); writes nothing
#   sbatch validate_pass1.sh regression \
#       $KRICO_RUNS/KRICO_0001/1993_11_15.nc \
#       $KRICO_RUNS/KRICO_0003/1994_01_15.nc
#
#   # the same with another ice source, against recruitment/data_nsidc/ or
#   # recruitment/data_glorys12_25km/
#   sbatch validate_pass1.sh regression_nsidc $KRICO_RUNS/KRICO_0001/1993_11_15.nc
#   sbatch validate_pass1.sh regression_glorys12_25km $KRICO_RUNS/KRICO_0001/1993_11_15.nc
#
#   # one-cohort test: classify one cohort into a test folder of
#   # recruitment/, then report on the file (outcome recomputed from the
#   # diagnostics, cross-checks, bytes per variable). A fourth word names
#   # another ice source.
#   sbatch validate_pass1.sh cohort $KRICO_RUNS/KRICO_0001/1993_11_15.nc data_v3.0_test
#   sbatch validate_pass1.sh cohort $KRICO_RUNS/KRICO_0001/1993_11_15.nc data_v3.0_test nsidc
#
#   # report on files already written
#   sbatch validate_pass1.sh report ../data_v3.0_glorys12/1993_11_15.nc
#
# Requires KRICO_GLORYS12_RELEASE_CELLS; with KRICO_GLORYS12 set, the
# regression also compares the sea ice at spawning with the monthly files.
# The other ice sources need KRICO_NSIDC or KRICO_GLORYS12_25KM.
# The report is the job output, logs/p1_val_<job id>.out; logs/ must exist.
# ============================================================================

set -euo pipefail

module load python3

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
RECRUITMENT_DIR="$(dirname "$SCRIPT_DIR")"

TEST="${1:?give a test: regression FILE..., cohort FILE FOLDER [SOURCE], or report FILE...}"
shift
case "${TEST}" in
    regression|report)
        [[ $# -gt 0 ]] || { echo "ERROR: ${TEST} needs files" >&2; exit 1; }
        python3 -u "${SCRIPT_DIR}/validate_pass1.py" "${TEST}" "$@"
        ;;
    regression_nsidc|regression_glorys12_25km)
        [[ $# -gt 0 ]] || { echo "ERROR: ${TEST} needs files" >&2; exit 1; }
        python3 -u "${SCRIPT_DIR}/validate_pass1.py" regression "$@" \
            --ice-source "${TEST#regression_}"
        ;;
    cohort)
        INPUT="${1:?cohort needs a trajectory file}"
        FOLDER="${2:?cohort needs a folder name under recruitment/, e.g. data_v3.0_test}"
        SOURCE="${3:-glorys12}"
        OUT_DIR="${RECRUITMENT_DIR}/${FOLDER}"
        mkdir -p "${OUT_DIR}"
        OUTPUT="${OUT_DIR}/$(basename "${INPUT}")"
        rm -f "${OUTPUT}"
        python3 -u "${SCRIPT_DIR}/process_cohort.py" "${INPUT}" "${OUTPUT}" \
            --ice-source "${SOURCE}"
        python3 -u "${SCRIPT_DIR}/validate_pass1.py" report "${OUTPUT}"
        ;;
    *)
        echo "ERROR: unknown test '${TEST}' (regression, regression_nsidc," \
             "regression_glorys12_25km, cohort or report)" >&2
        exit 1
        ;;
esac
