#!/bin/bash
#SBATCH --job-name=p2_check
#SBATCH --qos=nf
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/p2_check_%j.out
#SBATCH --error=logs/p2_check_%j.err

# ============================================================================
# Inputs of pass 2, checked for every cohort — one SLURM job
# ============================================================================
#
# Runs check_inputs_chl.py on the 32 spawning years (or on the years given):
# every cohort goes through the checks process_cohort_chl.py itself applies,
# and every problem is listed. Nothing is written. Only headers, time axes
# and release positions are read, so one job does all the years.
#
# Submit from recruitment/scripts, before run_year_chl.sh:
#   mkdir -p logs
#   sbatch check_inputs_chl.sh
#   sbatch check_inputs_chl.sh 1994          # one spawning year
# Then, with the job number sbatch printed:
#   grep -v ECMWF-INFO logs/p2_check_<job id>.out | grep "verdict\|PROBLEM\|WARNING\|Check\|CHECK"
#
# Requires KRICO_RUNS, KRICO_GLOBCOLOUR and KRICO_PISCES.
# ============================================================================

set -euo pipefail

module load python3

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

python3 -u "${SCRIPT_DIR}/check_inputs_chl.py" "$@"
