from pathlib import Path
import argparse

FILTER_LIBRARY = {
    "Filter_1": """\
*********** start of CM FLATFILT with identifier FILTER  ***********
0.75, RMAX
FILTER
8.1, ZMIN
1, NUMBER OF LAYERS
1, 0.1, # CONES, ZTHICK OF LAYER 1
0.75, 
0.75, 
0.521, 0.001, 0, 8, 
AL521ICRU
0.521, 0.001, 0, 9, 
AL521ICRU
*********** start of CM BLOCK with identifier COLLIMAT  ***********
""",

    "Filter_2": """\
*********** start of CM FLATFILT with identifier FILTER  ***********
0.75, RMAX
FILTER
8.1, ZMIN
1, NUMBER OF LAYERS
1, 0.2, # CONES, ZTHICK OF LAYER 1
0.75, 
0.75, 
0.521, 0.001, 0, 8, 
AL521ICRU
0.521, 0.001, 0, 9, 
AL521ICRU
*********** start of CM BLOCK with identifier COLLIMAT  ***********
""",

    "Filter_3": """\
*********** start of CM FLATFILT with identifier FILTER  ***********
0.75, RMAX
FILTER
8.1, ZMIN
1, NUMBER OF LAYERS
1, 0.25, # CONES, ZTHICK OF LAYER 1
0.75, 
0.75, 
0.521, 0.001, 0, 8, 
AL521ICRU
0.521, 0.001, 0, 9, 
AL521ICRU
*********** start of CM BLOCK with identifier COLLIMAT  ***********
""",

    "Filter_4": """\
*********** start of CM FLATFILT with identifier FILTER  ***********
0.75, RMAX
FILTER
8.1, ZMIN
2, NUMBER OF LAYERS
1, 0.05, # CONES, ZTHICK OF LAYER 1
0.75, 
0.75, 
1, 0.01, # CONES, ZTHICK OF LAYER 2
0.75, 
0.75, 
0.521, 0.001, 0, 8, 
AL521ICRU
0.521, 0.001, 0, 9, 
AL521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
*********** start of CM BLOCK with identifier COLLIMAT  ***********
""",

    "Filter_5": """\
*********** start of CM FLATFILT with identifier FILTER  ***********
0.75, RMAX
FILTER
8.1, ZMIN
2, NUMBER OF LAYERS
1, 0.1, # CONES, ZTHICK OF LAYER 1
0.75, 
0.75, 
1, 0.01, # CONES, ZTHICK OF LAYER 2
0.75, 
0.75, 
0.521, 0.001, 0, 8, 
AL521ICRU
0.521, 0.001, 0, 9, 
AL521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
*********** start of CM BLOCK with identifier COLLIMAT  ***********
""",

    "Filter_6": """\
*********** start of CM FLATFILT with identifier FILTER  ***********
0.75, RMAX
FILTER
8.1, ZMIN
2, NUMBER OF LAYERS
1, 0.15, # CONES, ZTHICK OF LAYER 1
0.75, 
0.75, 
1, 0.015, # CONES, ZTHICK OF LAYER 2
0.75, 
0.75, 
0.521, 0.001, 0, 8, 
AL521ICRU
0.521, 0.001, 0, 9, 
AL521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
*********** start of CM BLOCK with identifier COLLIMAT  ***********
""",

    "Filter_7": """\
*********** start of CM FLATFILT with identifier FILTER  ***********
0.75, RMAX
FILTER
8.1, ZMIN
2, NUMBER OF LAYERS
1, 0.1, # CONES, ZTHICK OF LAYER 1
0.75, 
0.75, 
1, 0.045, # CONES, ZTHICK OF LAYER 2
0.75, 
0.75, 
0.521, 0.001, 0, 8, 
AL521ICRU
0.521, 0.001, 0, 9, 
AL521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
*********** start of CM BLOCK with identifier COLLIMAT  ***********
""",

    "Filter_8": """\
*********** start of CM FLATFILT with identifier FILTER  ***********
0.75, RMAX
FILTER
8.1, ZMIN
2, NUMBER OF LAYERS
1, 0.1, # CONES, ZTHICK OF LAYER 1
0.75, 
0.75, 
1, 0.11, # CONES, ZTHICK OF LAYER 2
0.75, 
0.75, 
0.521, 0.001, 0, 8, 
AL521ICRU
0.521, 0.001, 0, 9, 
AL521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
*********** start of CM BLOCK with identifier COLLIMAT  ***********
""",

    "Filter_9": """\
*********** start of CM FLATFILT with identifier FILTER  ***********
0.75, RMAX
FILTER
8.1, ZMIN
3, NUMBER OF LAYERS
1, 0.15, # CONES, ZTHICK OF LAYER 1
0.75, 
0.75, 
1, 0.025, # CONES, ZTHICK OF LAYER 2
0.75, 
0.75, 
1, 0.05, # CONES, ZTHICK OF LAYER 3
0.75, 
0.75, 
0.521, 0.001, 0, 8, 
AL521ICRU
0.521, 0.001, 0, 9, 
AL521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
0.521, 0.001, 0, 0, 
CU521ICRU
0.521, 0.001, 0, 0, 
SN700ICRU
0.521, 0.001, 0, 0, 
SN700ICRU
*********** start of CM BLOCK with identifier COLLIMAT  ***********
"""
}

# --------------------------------------------------
# Core functions
# --------------------------------------------------
def swap_mono_energy(lines, energy_MeV):
    out = []
    replace_next = False
    replaced = False

    for line in lines:
        if replace_next:
            out.append(f"{energy_MeV:.3f}\n")
            replace_next = False
            replaced = True
            continue

        if line.strip() == "0, MONOENERGETIC":
            replace_next = True
            out.append(line)
            continue

        out.append(line)

    if not replaced:
        raise RuntimeError("Monoenergetic energy not found.")

    return out


def swap_filter_block(lines, filter_block_text):
    out = []
    inside_filter = False
    replaced = False

    for line in lines:
        if "*********** start of CM FLATFILT with identifier FILTER  ***********" in line:
            inside_filter = True
            out.append(filter_block_text)
            replaced = True
            continue

        if inside_filter:
            if "*********** start of CM BLOCK with identifier COLLIMAT  ***********" in line:
                inside_filter = False
            continue

        out.append(line)

    if not replaced:
        raise RuntimeError("FILTER block not found.")

    return out


# --------------------------------------------------
# CLI
# --------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Modify BEAMnrc input files (energy, filter)."
    )

    parser.add_argument(
        "-e", "--energy",
        type=float,
        required=True,
        help="Monoenergetic energy in MeV (e.g. 0.080)"
    )

    parser.add_argument(
        "-f", "--filter",
        required=True,
        choices=FILTER_LIBRARY.keys(),
        help="Filter name to swap in"
    )

    args = parser.parse_args()

    # Master Template Path
    templatepath = Path("template_files/template.egsinp")

    lines = templatepath.read_text().splitlines(keepends=True)

    if args.energy is not None:
        lines = swap_mono_energy(lines, args.energy)

    if args.filter is not None:
        lines = swap_filter_block(lines, FILTER_LIBRARY[args.filter])

    output = Path(f"{args.energy*1000:.0f}kvP_{args.filter}.egsinp")
    
    output.write_text("".join(lines))
    print(f"Output file {output} written.")

if __name__ == "__main__":
    main()
