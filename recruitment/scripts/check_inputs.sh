#!/bin/bash
#SBATCH --job-name=p1_check
#SBATCH --qos=nf
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --array=1-32
#SBATCH --output=logs/p1_check_%A_%a.out
#SBATCH --error=logs/p1_check_%A_%a.err

# ============================================================================
# Inputs of pass 1, checked for every cohort — SLURM array
# ============================================================================
#
# One array task per spawning year, index 1-32 = 1994-2025, as in
# run_year.sh. Each task runs check_inputs.py on its year: every cohort goes
# through the checks process_cohort.py itself applies, and every problem is
# listed. Nothing is written.
#
# Submit from recruitment/scripts, before run_year.sh:
#   mkdir -p logs
#   sbatch check_inputs.sh
#   sbatch check_inputs.sh nsidc             # before a run with NSIDC sea ice
#   sbatch check_inputs.sh glorys12_25km     # or GLORYS12 ice on the NSIDC grid
# Then, with the job number sbatch printed:
#   grep -h "verdict\|PROBLEM\|WARNING" logs/p1_check_<job id>_*.out
#
# Requires KRICO_RUNS and KRICO_GLORYS12_RELEASE_CELLS, and KRICO_NSIDC or
# KRICO_GLORYS12_25KM for the other ice sources.
# ============================================================================

set -euo pipefail

module load python3

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
INDEX=${SLURM_ARRAY_TASK_ID:?run through sbatch, or call check_inputs.py YEAR directly}
YEAR=$((1993 + INDEX))

python3 -u "${SCRIPT_DIR}/check_inputs.py" "${YEAR}" "$@"
