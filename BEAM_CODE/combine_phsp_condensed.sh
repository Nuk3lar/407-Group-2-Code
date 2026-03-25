#!/usr/bin/env sh   
set -e
if [ "$#" -ne 2 ]; then
    echo "Usage: combine_phsp.sh <input> <nproc>"
    exit 1
fi
INPUT_FILE="$1"
NPROC="$2"
add_phsp="$HEN_HOUSE/bin/linux64/addphsp"
# Load known-good environment
if [ ! -f "./egsnrc_env.txt" ]; then
    echo "ERROR: egsnrc_env.txt not found. Run 'env > egsnrc_env.txt' first."
    exit 2
fi
# Export environment
while IFS='=' read -r key value; do
    export "$key=$value"
done < ./egsnrc_env.txt
# Run
$add_phsp $INPUT_FILE "$INPUT_FILE"_combined $NPROC 1 1 0
$add_phsp $INPUT_FILE "$INPUT_FILE"_combined $NPROC 1 2 0
$add_phsp $INPUT_FILE "$INPUT_FILE"_combined $NPROC 1 3 0
