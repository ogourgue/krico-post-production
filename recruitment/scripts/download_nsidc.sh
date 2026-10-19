#!/bin/bash
#SBATCH --job-name=nsidc_dl
#SBATCH --qos=nf
#SBATCH --time=02:00:00
#SBATCH --mem=2G
#SBATCH --array=0-32%6
#SBATCH --output=logs/nsidc_dl_%A_%a.out
#SBATCH --error=logs/nsidc_dl_%A_%a.err

# ============================================================================
# Download the daily NSIDC G02202 v6 sea-ice files — SLURM array
# ============================================================================
#
# One array task per calendar year: index 0-32 maps to 1993-2025 (SLURM caps
# array indices below MaxArraySize, so the year itself cannot be the index).
# At most 6 tasks run at once so that NSIDC does not see 33 parallel clients.
# Each task calls download_nsidc.py for its year; the script is resumable, so
# a failed task can simply be resubmitted with index = year - 1993:
#   sbatch --array=14 download_nsidc.sh      # 2007
#
# When all tasks have finished, write the provenance file (checksums and the
# dates for which the archive offers no file), interactively, it takes under a
# minute:
#   python3 download_nsidc.py --provenance
# and commit nsidc_daily_provenance.txt alongside this script.
#
# Requires KRICO_NSIDC (e.g. /scratch/cvan/KRICO/Pre/NSIDC). Job logs go to
# recruitment/scripts/logs/, which must exist at submission time:
#   mkdir -p logs
#   sbatch download_nsidc.sh
# ============================================================================

set -euo pipefail

module load python3

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
INDEX=${SLURM_ARRAY_TASK_ID:?run through sbatch, or call download_nsidc.py YEAR directly}
YEAR=$((1993 + INDEX))

python3 "${SCRIPT_DIR}/download_nsidc.py" "${YEAR}"
