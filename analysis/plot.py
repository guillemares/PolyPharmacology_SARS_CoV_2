from 3D_RNAlysis import rna


def plot_ensemble_interactions(rnas):
    """
    given an ensemble of RNA structures (all same bases) coming from MD, ML, ...
    """


    return 0


if if __name__ == "__main__":
    import glob
    rnas_files = glob.glob('*.pdb')
    rnas = []
    for rna_file in rnas_files:
        rnaobj = rna.RNA(pdbfile=rna_file)
        rnas.append(rnaobj)
