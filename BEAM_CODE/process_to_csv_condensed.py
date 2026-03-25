#!/usr/bin/env python3

import sys
from pathlib import Path

HEADER_LINES = 17   # number of lines to strip (fixed)

def main():
    if len(sys.argv) != 2:
        print("Usage: python convert_to_csv.py <input_file>")
        sys.exit(1)

    input_file = Path(sys.argv[1])

    if not input_file.exists():
        print(f"Error: file not found: {input_file}")
        sys.exit(1)

    output_file_csv = input_file.with_suffix(".csv")

    # Read input file
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Remove header
    data_lines = lines[HEADER_LINES:]

    # Write CSV
    with open(output_file_csv, "w") as fout:
        for line in data_lines:
            if line.strip():            # skip blank lines
                cols = line.split()     # split on whitespace
                fout.write(",".join(cols) + "\n")

    print(f"CSV written to: {output_file_csv}")

if __name__ == "__main__":
    main()
