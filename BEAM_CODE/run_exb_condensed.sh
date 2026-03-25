#!/usr/bin/env sh   
set -e
if [ "$#" -ne 4 ]; then
    echo "Usage: run_exb.sh <beam_code> <input.egsinp> <pegs.pegs4dat> <nproc>"
    exit 1
fi
BEAM_CODE="$1"
INPUT_FILE="$2"
PEGS_FILE="$3"
NPROC="$4"
EXB="$HEN_HOUSE/scripts/run_user_code_batch"
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
exec $EXB "$BEAM_CODE" "$INPUT_FILE" "$PEGS_FILE" "p=$NPROC"
