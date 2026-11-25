#!/bin/bash
#SBATCH --job-name=cells_pre
#SBATCH --qos=nf
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --array=1-32
#SBATCH --output=logs/cells_pre_%A_%a.out
#SBATCH --error=logs/cells_pre_%A_%a.err

# ============================================================================
# GLORYS12 sea ice and temperature at the release cells — SLURM array
# ============================================================================
#
# One array task per spawning year, index 1-32 = 1994-2025, as in
# run_year.sh. Each task writes glorys12_release_cells_YYYY.nc into
# $KRICO_GLORYS12_RELEASE_CELLS from the monthly sea-ice and temperature
# files in $KRICO_GLORYS12 (22 October of the previous year to 30 September).
# A year whose file exists is skipped, so a failed task is resubmitted with
# index = year - 1993:
#   sbatch --array=14 preprocess_glorys12_release_cells.sh    # 2007
#
# The files record the commit that produced them: the task stops if tracked
# files are modified or if the commit is not pushed.
#
# Check the inputs once before submitting (login node, a few minutes):
#   python3 preprocess_glorys12_release_cells.py --check
#
# Requires KRICO_GLORYS12 and KRICO_GLORYS12_RELEASE_CELLS. Logs go to
# recruitment/scripts/logs/, which must exist:
#   mkdir -p logs
#   sbatch --array=1 preprocess_glorys12_release_cells.sh     # 1994, as a test
#   sbatch --array=2-32 preprocess_glorys12_release_cells.sh  # the rest
# ============================================================================

set -euo pipefail

module load python3

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
INDEX=${SLURM_ARRAY_TASK_ID:?run through sbatch, or call preprocess_glorys12_release_cells.py YEAR directly}
YEAR=$((1993 + INDEX))

python3 -u "${SCRIPT_DIR}/preprocess_glorys12_release_cells.py" "${YEAR}"
