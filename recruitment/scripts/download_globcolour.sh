#!/bin/bash
#SBATCH --job-name=globcolour_dl
#SBATCH --qos=nf
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/globcolour_dl_%j.out
#SBATCH --error=logs/globcolour_dl_%j.err

# ============================================================================
# Download the GlobColour daily chlorophyll climatology — single SLURM job
# ============================================================================
#
# One file of 366 days, so one job and no array. The script does nothing when
# the file is already complete, so a failed job can simply be resubmitted.
#
# When the job has finished, write the provenance file, interactively:
#   python3 download_globcolour.py --provenance
# and commit globcolour_chl_provenance.txt alongside this script.
#
# Requires:
#   KRICO_GLOBCOLOUR   e.g. /scratch/cvan/KRICO/Pre/GLOBCOLOUR
#   a Copernicus Marine login stored on the machine (copernicusmarine login)
# The python3 module of the ECMWF HPC provides the Copernicus Marine toolbox.
# Where it does not, set KRICO_CMT_VENV to a virtual environment that holds
# it (e.g. $PERM/venvs/copernicusmarine) and it is activated here.
# Job logs go to recruitment/scripts/logs/, which must exist at submission:
#   mkdir -p logs
#   sbatch download_globcolour.sh
# ============================================================================

set -euo pipefail

module load python3
if [[ -n "${KRICO_CMT_VENV:-}" ]]; then
    # The activation script is not written for "set -u".
    set +u
    source "${KRICO_CMT_VENV}/bin/activate"
    set -u
fi

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# -u: unbuffered, so that the log is written as the job runs.
python3 -u "${SCRIPT_DIR}/download_globcolour.py"
