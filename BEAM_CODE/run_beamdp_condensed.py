import os
import subprocess

# -------------------------
# USER SETTINGS
# -------------------------

# Energies in MeV
energies = [0.080, 0.100, 0.120, 0.150, 0.180, 0.200, 0.250, 0.300]

emin = 0.001
nbins = 2000
output_dir = "/home/meepn/EGSnrc/EGS_HOME/BEAM_XStrahl300_2/spectra2"
input_dir  = "/home/meepn/EGSnrc/EGS_HOME/BEAM_XStrahl300_2/beamdp_inputs2"

run_beamdp = True   # set True if you want auto-run

# -------------------------
# SETUP
# -------------------------

os.makedirs(output_dir, exist_ok=True)
os.makedirs(input_dir, exist_ok=True)

# -------------------------
# LOOP OVER ENERGIES
# -------------------------

for emax in energies:
    i = 2
    while i < 3:
        z = 1
        while z < 2:
            tag = f"{int(emax*1000)}kvP_Filter_{z}"

            input_file = os.path.join(input_dir, f"beamdp_{tag}_{i}.input")
            output_file = os.path.join(output_dir, f"{tag}_{i}.txt")
            output_file_csv = os.path.join(output_dir, f"{tag}_{i}.csv")
            
            phsp_file = f"/home/meepn/EGSnrc/EGS_HOME/BEAM_XStrahl300_2/{tag}_combined.egsphsp{i}"
            beamdp_input = f"""BeamDP Input File
3
0, 0, 10
{nbins}, {emin}, {emax}
0, 0, 0
{phsp_file}
{output_file}
1
1
0
0
"""

            with open(input_file, "w") as f:
                f.write(beamdp_input)

            print(f"Created: {input_file}")

            if run_beamdp:
                subprocess.run(
                    ["beamdp"],
                    stdin=open(input_file),
                    cwd=os.path.dirname(phsp_file)
                )
                with open(output_file, "r") as f:
                    lines = f.readlines()

                with open(output_file, "w") as f:
                    f.writelines(lines[17:])
                
                with open(output_file, "r") as fin, open(output_file_csv, "w") as fout:
                    for line in fin:
                        if line.strip():  # skip blank lines
                            cols = line.split()   # split on any whitespace
                            fout.write(",".join(cols) + "\n")
            z = z + 1
                
        i = i + 1
