#!/bin/bash
#
# Archive recruitment data by spawning year
# ==========================================
#
# Groups recruitment files (one per release date) into spawning-year archives.
# Each spawning year includes releases from Nov 15 of the previous calendar year
# through Mar 14 of the spawning year (i.e., the four release months: Nov, Dec, Jan, Feb).
#
# Spawning year 1994 → files 1993_11_15.nc through 1994_03_14.nc
# Spawning year 1995 → files 1994_11_15.nc through 1995_03_14.nc
# ...
# Spawning year 2025 → files 2024_11_15.nc through 2025_03_14.nc
#
# Assumes this script is in recruitment/scripts/ and creates archives in recruitment/archives/
#
# Usage
# -----
#   cd recruitment/scripts
#   ./archive_by_year.sh
#
# Output
# ------
#   ../archives/recruitment_1994.tar.gz
#   ../archives/recruitment_1995.tar.gz
#   ...
#   ../archives/recruitment_2025.tar.gz
#

set -euo pipefail

# Fixed directory structure: script is in recruitment/scripts/
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RECRUITMENT_DIR="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="$RECRUITMENT_DIR/data"
OUTPUT_DIR="$RECRUITMENT_DIR/archives"

# Validate input directory
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: recruitment data directory not found: $INPUT_DIR" >&2
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "================================================================"
echo "Archiving recruitment files by spawning year"
echo "================================================================"
echo "Input directory:  $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo

# Count files for sanity check
file_count=$(find "$INPUT_DIR" -maxdepth 1 -name "[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9].nc" | wc -l)
echo "Found $file_count recruitment files"
echo

# Loop over spawning years (1994-2025)
for spawning_year in {1994..2025}; do
    prev_year=$((spawning_year - 1))

    echo -n "Spawning year $spawning_year (files from ${prev_year}_11 through ${spawning_year}_03)... "

    # Build list of files for this spawning year
    # Includes:
    #   - ${prev_year}_11_*.nc and ${prev_year}_12_*.nc (Nov and Dec of previous year)
    #   - ${spawning_year}_01_*.nc, ${spawning_year}_02_*.nc, ${spawning_year}_03_*.nc (Jan, Feb, Mar of spawning year)

    files_to_archive=()

    # November of previous year
    while IFS= read -r -d '' file; do
        files_to_archive+=("$file")
    done < <(find "$INPUT_DIR" -maxdepth 1 -name "${prev_year}_11_[0-9][0-9].nc" -print0)

    # December of previous year
    while IFS= read -r -d '' file; do
        files_to_archive+=("$file")
    done < <(find "$INPUT_DIR" -maxdepth 1 -name "${prev_year}_12_[0-9][0-9].nc" -print0)

    # January of spawning year
    while IFS= read -r -d '' file; do
        files_to_archive+=("$file")
    done < <(find "$INPUT_DIR" -maxdepth 1 -name "${spawning_year}_01_[0-9][0-9].nc" -print0)

    # February of spawning year
    while IFS= read -r -d '' file; do
        files_to_archive+=("$file")
    done < <(find "$INPUT_DIR" -maxdepth 1 -name "${spawning_year}_02_[0-9][0-9].nc" -print0)

    # March of spawning year (all files; last release is typically 03_14 by design)
    while IFS= read -r -d '' file; do
        files_to_archive+=("$file")
    done < <(find "$INPUT_DIR" -maxdepth 1 -name "${spawning_year}_03_[0-9][0-9].nc" -print0)

    # Sort files chronologically
    IFS=$'\n' files_to_archive=($(sort <<<"${files_to_archive[*]}"))
    unset IFS

    if [[ ${#files_to_archive[@]} -eq 0 ]]; then
        echo "SKIP (no files found)"
        continue
    fi

    # Create tar.gz archive
    archive_name="recruitment_${spawning_year}.tar.gz"
    archive_path="$OUTPUT_DIR/$archive_name"

    # Use relative paths in archive for cleaner extraction
    # Change to input directory, archive files by passing basenames via stdin
    cd "$INPUT_DIR"
    {
        for file in "${files_to_archive[@]}"; do
            basename "$file"
        done
    } | tar --gzip --create --file "$archive_path" -T -
    cd - > /dev/null

    # Report size and count
    archive_size=$(du -h "$archive_path" | cut -f1)
    file_count=${#files_to_archive[@]}

    echo "OK ($file_count files, $archive_size)"
done

echo
echo "================================================================"
echo "Archiving complete"
echo "================================================================"
echo
echo "Archives created in: $OUTPUT_DIR"
echo
ls -lh "$OUTPUT_DIR"/recruitment_*.tar.gz | awk '{print "  " $9 " (" $5 ")"}'
echo
echo "Total size:"
du -sh "$OUTPUT_DIR" | awk '{print "  " $1}'
echo

echo "Verification: list contents of first archive"
echo "  tar -tzf $OUTPUT_DIR/recruitment_1994.tar.gz | head -5"
tar -tzf "$OUTPUT_DIR/recruitment_1994.tar.gz" | head -5
echo "  ..."