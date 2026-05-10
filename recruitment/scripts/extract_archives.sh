#!/bin/bash
#
# Extract recruitment archives into data folder
# =============================================
#
# Extracts all spawning-year tar.gz archives from recruitment/archives/
# into recruitment/data/. Existing files are skipped.
#
# Assumes this script is in recruitment/scripts/
#
# Usage
# -----
#   cd recruitment/scripts
#   ./extract_archives.sh
#

set -euo pipefail

# Fixed directory structure: script is in recruitment/scripts/
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RECRUITMENT_DIR="$(dirname "$SCRIPT_DIR")"
ARCHIVES_DIR="$RECRUITMENT_DIR/archives"
DATA_DIR="$RECRUITMENT_DIR/data"

# Validate archives directory
if [[ ! -d "$ARCHIVES_DIR" ]]; then
    echo "ERROR: archives directory not found: $ARCHIVES_DIR" >&2
    exit 1
fi

# Create data directory if needed
mkdir -p "$DATA_DIR"

echo "================================================================"
echo "Extracting recruitment archives"
echo "================================================================"
echo "Archives: $ARCHIVES_DIR"
echo "Data:     $DATA_DIR"
echo

# Count available archives
archive_count=$(find "$ARCHIVES_DIR" -maxdepth 1 -name "recruitment_[0-9][0-9][0-9][0-9].tar.gz" | wc -l)
echo "Found $archive_count archives"
echo

n_extracted=0
n_skipped=0

for archive in $(find "$ARCHIVES_DIR" -maxdepth 1 -name "recruitment_[0-9][0-9][0-9][0-9].tar.gz" | sort); do
    archive_name=$(basename "$archive")
    year="${archive_name:12:4}"

    echo -n "Spawning year $year... "

    # Check if files from this archive already exist in data/
    first_file=$(tar -tzf "$archive" | head -1)
    if [[ -f "$DATA_DIR/$first_file" ]]; then
        echo "SKIP (already extracted)"
        n_skipped=$((n_skipped + 1))
        continue
    fi

    # Extract into data directory
    tar --gzip --extract --file "$archive" --directory "$DATA_DIR"
    file_count=$(tar -tzf "$archive" | wc -l)
    echo "OK ($file_count files)"
    n_extracted=$((n_extracted + 1))
done

echo
echo "================================================================"
echo "Done"
echo "  extracted: $n_extracted archives"
echo "  skipped:   $n_skipped archives (already extracted)"
echo "================================================================"