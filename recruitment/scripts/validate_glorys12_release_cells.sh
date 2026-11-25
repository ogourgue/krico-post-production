#!/bin/bash
#SBATCH --job-name=cells_val
#SBATCH --qos=nf
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/cells_val_%j.out
#SBATCH --error=logs/cells_val_%j.err

# ============================================================================
# Acceptance tests of the release-cell extraction — SLURM driver
# ============================================================================
#
# Runs validate_glorys12_release_cells.py with the arguments given after the
# script name, off the login node. Submit from recruitment/scripts:
#
#   # coverage, over the classification files of the whole dataset
#   sbatch validate_glorys12_release_cells.sh --releases ../data
#
#   # day 0, once the extraction of the cohorts' spawning year exists
#   sbatch validate_glorys12_release_cells.sh \
#       $KRICO_RUNS/KRICO_0001/1993_11_15.nc \
#       $KRICO_RUNS/KRICO_0003/1994_01_15.nc \
#       $KRICO_RUNS/KRICO_0004/1994_03_14.nc
#
# Requires KRICO_GLORYS12, and KRICO_GLORYS12_RELEASE_CELLS for day 0. The
# report is the job output, logs/cells_val_<job id>.out; logs/ must exist.
# ============================================================================

set -euo pipefail

module load python3

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

python3 -u "${SCRIPT_DIR}/validate_glorys12_release_cells.py" "$@"
