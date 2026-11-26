#!/bin/bash
#SBATCH --job-name=p2_check
#SBATCH --qos=nf
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --array=1-32
#SBATCH --output=logs/p2_check_%A_%a.out
#SBATCH --error=logs/p2_check_%A_%a.err

# ============================================================================
# Inputs of pass 2, checked for every cohort — SLURM array
# ============================================================================
#
# One array task per spawning year, index 1-32 = 1994-2025, as in
# run_year_chl.sh. Each task runs check_inputs_chl.py on its year: every
# cohort goes through the checks process_cohort_chl.py itself applies, and
# every problem is listed. Nothing is written.
#
# Submit from recruitment/scripts, before run_year_chl.sh:
#   mkdir -p logs
#   sbatch check_inputs_chl.sh
# Then, with the job number sbatch printed:
#   grep -h "verdict\|PROBLEM\|WARNING" logs/p2_check_<job id>_*.out
#
# Requires KRICO_RUNS, KRICO_GLOBCOLOUR and KRICO_PISCES.
# ============================================================================

set -euo pipefail

module load python3

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

INDEX=${SLURM_ARRAY_TASK_ID:?run through sbatch, or call check_inputs_chl.py YEAR directly}
YEAR=$((1993 + INDEX))

python3 -u "${SCRIPT_DIR}/check_inputs_chl.py" "${YEAR}"
