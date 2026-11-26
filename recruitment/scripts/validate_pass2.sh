#!/bin/bash
#SBATCH --job-name=p2_val
#SBATCH --qos=nf
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/p2_val_%j.out
#SBATCH --error=logs/p2_val_%j.err

# ============================================================================
# Tests of pass 2 on real cohorts — SLURM driver
# ============================================================================
#
# The test is named by a plain word, because ECMWF's sbatch refuses script
# arguments that start with a dash. Submit from recruitment/scripts:
#
#   # one-cohort test: compute the chlorophyll indices of one cohort into a
#   # test folder of recruitment/, then report on the file (consistency,
#   # windows against the pass 1 file in recruitment/data_v3.0_glorys12/,
#   # recomputation of a sample by another route, bytes per variable)
#   sbatch validate_pass2.sh cohort $KRICO_RUNS/KRICO_0003/1994_01_15.nc data_v3.0_chl_test
#
#   # report on files already written
#   sbatch validate_pass2.sh report ../data_v3.0_chl/1994_01_15.nc
#
# Requires KRICO_RUNS, KRICO_GLOBCOLOUR and KRICO_PISCES.
# The report is the job output, logs/p2_val_<job id>.out; logs/ must exist.
# ============================================================================

set -euo pipefail

module load python3

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
RECRUITMENT_DIR="$(dirname "$SCRIPT_DIR")"

TEST="${1:?give a test: cohort FILE FOLDER, or report FILE...}"
shift
case "${TEST}" in
    report)
        [[ $# -gt 0 ]] || { echo "ERROR: report needs files" >&2; exit 1; }
        python3 -u "${SCRIPT_DIR}/validate_pass2.py" "$@"
        ;;
    cohort)
        INPUT="${1:?cohort needs a trajectory file}"
        FOLDER="${2:?cohort needs a folder name under recruitment/, e.g. data_v3.0_chl_test}"
        OUT_DIR="${RECRUITMENT_DIR}/${FOLDER}"
        mkdir -p "${OUT_DIR}"
        OUTPUT="${OUT_DIR}/$(basename "${INPUT}")"
        rm -f "${OUTPUT}"
        python3 -u "${SCRIPT_DIR}/process_cohort_chl.py" "${INPUT}" "${OUTPUT}"
        python3 -u "${SCRIPT_DIR}/validate_pass2.py" "${OUTPUT}"
        ;;
    *)
        echo "ERROR: unknown test '${TEST}' (cohort or report)" >&2
        exit 1
        ;;
esac
