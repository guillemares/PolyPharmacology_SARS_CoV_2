#!/usr/bin/env python
import os
import glob

"""
In this script we split pdb files by models, since some
of the 3D structure prediction tools generate multiple
models in a single pdb file and the tools in aRNAlysis
pipeline require single model pdb files.
"""

indir = "path/to/pdb_files" # Change this to the path
                            # containing the pdb files.


def split_pdbs(file, outdir):
    with open(file, "r") as f:
        lines = f.readlines()
        model = 1
        for i in range(len(lines)):
            if "MODEL" in lines[i]:
                output_file = os.path.join(outdir, file.split(".pdb")[0] + "_MODEL" + str(model) + ".pdb")
                with open(output_file, "w") as g:
                    g.write(lines[i])
                    j = i + 1
                    while "ENDMDL" not in lines[j]:
                        g.write(lines[j])
                        j += 1
                    g.write(lines[j])
                    model += 1
        with open(file, "a") as h:
            h.write("END\n")

    return 0

def main():
    outdir = os.path.join(indir, "Models")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    pdb_files = glob.glob(os.path.join(indir, "*.pdb"))
    for pdb_file in pdb_files:
        split_pdbs(pdb_file, outdir)

    return 0
