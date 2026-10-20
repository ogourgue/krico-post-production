#!/bin/bash
#SBATCH --job-name=nsidc_pre
#SBATCH --qos=nf
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --array=0-33
#SBATCH --output=logs/nsidc_pre_%A_%a.out
#SBATCH --error=logs/nsidc_pre_%A_%a.err

# ============================================================================
# Stack the daily NSIDC G02202 v6 files into monthly files — SLURM array
# ============================================================================
#
# One array task per calendar year, index 0-33 = 1993-2026, as in
# download_nsidc.sh. Each task writes the nsidc_ice_YYYY_MM.nc files of its
# year into $KRICO_NSIDC; months already there are skipped, so a failed task
# is resubmitted with index = year - 1993:
#   sbatch --array=14 preprocess_nsidc.sh    # 2007
#
# Run the grid check once before submitting (a few seconds, on the login node):
#   python3 preprocess_nsidc.py --check-grid
#
# Requires KRICO_NSIDC. Logs go to recruitment/scripts/logs/, which must exist:
#   mkdir -p logs
#   sbatch preprocess_nsidc.sh
# ============================================================================

set -euo pipefail

module load python3

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
INDEX=${SLURM_ARRAY_TASK_ID:?run through sbatch, or call preprocess_nsidc.py YEAR directly}
YEAR=$((1993 + INDEX))

python3 "${SCRIPT_DIR}/preprocess_nsidc.py" "${YEAR}"
