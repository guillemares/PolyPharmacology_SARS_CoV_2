from aRNAlysis import rna
import glob
import doctest
import argparse

def plot_ensemble_interactions(rnas):
    """
    given an ensemble of RNA structures (all same bases) coming from MD, ML, ...
    """


    return 0


if __name__ == "__main__":
    rnas_files = glob.glob('*.pdb')
    rnas = []
    for rna_file in rnas_files:
        rnaobj = rna.RNA(pdbfile=rna_file)
        rnas.append(rnaobj)
