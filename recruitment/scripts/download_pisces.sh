#!/bin/bash
#SBATCH --job-name=pisces_dl
#SBATCH --qos=nf
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --array=0-32%6
#SBATCH --output=logs/pisces_dl_%A_%a.out
#SBATCH --error=logs/pisces_dl_%A_%a.err

# ============================================================================
# Download the daily 3D chlorophyll of the PISCES hindcast — SLURM array
# ============================================================================
#
# One array task per calendar year: index 0-32 maps to 1993-2025, as in
# download_nsidc.sh. At most 6 tasks run at once so that Copernicus Marine
# does not see 33 parallel clients. Each task calls download_pisces.py for its
# year; the script is resumable, so a failed task can simply be resubmitted
# with index = year - 1993:
#   sbatch --array=14 download_pisces.sh      # 2007
#
# When all tasks have finished, write the provenance file (dataset version,
# toolbox version, days and size of every file), interactively:
#   python3 download_pisces.py --provenance
# and commit pisces_chl_provenance.txt alongside this script.
#
# Requires:
#   KRICO_PISCES   e.g. /scratch/cvan/KRICO/Pre/PISCES
#   KRICO_CMT_VENV virtual environment holding the Copernicus Marine toolbox,
#                  e.g. $PERM/venvs/copernicusmarine
#   a Copernicus Marine login stored on the machine (copernicusmarine login)
# Job logs go to recruitment/scripts/logs/, which must exist at submission:
#   mkdir -p logs
#   sbatch download_pisces.sh
# ============================================================================

set -euo pipefail

module load python3
# The activation script is not written for "set -u".
set +u
source "${KRICO_CMT_VENV:?set KRICO_CMT_VENV to the toolbox virtual environment}/bin/activate"
set -u

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
INDEX=${SLURM_ARRAY_TASK_ID:?run through sbatch, or call download_pisces.py YEAR directly}
YEAR=$((1993 + INDEX))

python3 "${SCRIPT_DIR}/download_pisces.py" "${YEAR}"
