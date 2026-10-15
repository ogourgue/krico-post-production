#!/bin/bash
#SBATCH --job-name=g25_pre
#SBATCH --qos=nf
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --array=0-32
#SBATCH --output=logs/g25_pre_%A_%a.out
#SBATCH --error=logs/g25_pre_%A_%a.err

# ============================================================================
# Bin GLORYS12 ice onto the NSIDC 25 km grid, monthly files — SLURM array
# ============================================================================
#
# One array task per calendar year, index 0-32 = 1993-2025, as in
# preprocess_nsidc.sh. Each task writes the glorys12_25km_ice_YYYY_MM.nc files
# of its year into $KRICO_GLORYS12_25KM from the GLORYS12 monthly ice in
# $KRICO_GLORYS12, masked with the NSIDC ancillary in $KRICO_NSIDC; months
# already there are skipped, so a failed task is resubmitted with
# index = year - 1993:
#   sbatch --array=14 preprocess_glorys12_25km.sh    # 2007
#
# Run the static-map check once before submitting (seconds, login node):
#   python3 preprocess_glorys12_25km.py --check
#
# Requires KRICO_GLORYS12, KRICO_NSIDC, KRICO_GLORYS12_25KM. Logs go to
# recruitment/scripts/logs/, which must exist:
#   mkdir -p logs
#   sbatch preprocess_glorys12_25km.sh
# ============================================================================

set -euo pipefail

module load python3

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
INDEX=${SLURM_ARRAY_TASK_ID:?run through sbatch, or call preprocess_glorys12_25km.py YEAR directly}
YEAR=$((1993 + INDEX))

python3 "${SCRIPT_DIR}/preprocess_glorys12_25km.py" "${YEAR}"
