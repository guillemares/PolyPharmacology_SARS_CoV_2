#!/usr/bin/env python
import os
import glob

def nt_terminal(pdb_file):
    r"""
    When using 3dRNA structures, the terminal nucleotides are not
    properly defined for the superimposition. This function will
    fix the terminal nucleotides of the 3dRNA pdb files.
    """

    with open(pdb_file, "r") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if (line.startswith("ATOM") or line.startswith("TER")) and (line[18:20] == "G5" or line[18:20] == "G3"):
            lines[i] = line[:18] + " G" + line[20:]
        elif (line.startswith("ATOM") or line.startswith("TER")) and (line[18:20] == "C5" or line[18:20] == "C3"):
            lines[i] = line[:18] + " C" + line[20:]
        elif (line.startswith("ATOM") or line.startswith("TER")) and (line[18:20] == "A5" or line[18:20] == "A3"):
            lines[i] = line[:18] + " A" + line[20:]
        elif (line.startswith("ATOM") or line.startswith("TER")) and (line[18:20] == "U5" or line[18:20] == "U3"):
            lines[i] = line[:18] + " U" + line[20:]

    with open(pdb_file, "w") as file:
        file.writelines(lines)
