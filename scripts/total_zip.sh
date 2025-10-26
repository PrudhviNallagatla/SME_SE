#!/bin/bash
#
# This script compresses the entire 'raw_output' directory into a single
# .tar.zst file for easy downloading.
#
# It uses 'zstd' for fast, multi-threaded compression.
#
# IMPORTANT: Run this *from* a VM with your /home/rimuru/workspace
#            file share already mounted!
#
set -euo pipefail

WORKSPACE_ROOT="/home/rimuru/workspace"
SOURCE_DIR="${WORKSPACE_ROOT}/data/raw_output"
OUTPUT_FILE="${WORKSPACE_ROOT}/zips/simulation_results_$(date +%Y-%m-%d).tar.zst"

echo "=================================================="
echo "Starting Compression..."
echo "=================================================="
echo "Source:  ${SOURCE_DIR}"
echo "Output:  ${OUTPUT_FILE}"
echo ""
echo "This may take a long time for 150GB+ of data..."

# Check if the source directory exists
if [ ! -d "${SOURCE_DIR}" ]; then
    echo "ERROR: Source directory not found!"
    echo "${SOURCE_DIR}"
    exit 1
fi

# Run the compression:
# 'tar' options:
#   -c : Create an archive
#   -f : Use the specified archive FILE
#   -I : Pipe through a compression Program
#
# 'zstd' options:
#   -T0 : Use all available CPU cores (threads)
#
# We change directory to 'data/' so the archive paths
# don't include the full '/home/rimuru/workspace'
cd "${WORKSPACE_ROOT}/data"
tar -cf "${OUTPUT_FILE}" -I 'zstd -T0' "raw_output"

# Go back to the original directory
cd "${WORKSPACE_ROOT}"

echo ""
echo "=================================================="
echo "Compression Finished!"
echo "=================================================="

# Final check
if [ -f "${OUTPUT_FILE}" ]; then
    echo "✓ Success! Archive created at:"
    echo "  ${OUTPUT_FILE}"
    ls -lh "${OUTPUT_FILE}"
else
    echo "✗ FAILED! Archive was not created."
    exit 1
fi