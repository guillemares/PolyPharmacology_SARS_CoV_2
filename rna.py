#!/usr/bin/env python
import biotite.structure.io.pdb as pdb
import biotite.structure as struc
from biotite.structure.error import BadStructureError
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os


class RNA(object):
    """
    """
    def __init__(self, pdbfile, init_index=None, test=False, outdir=None):
        """
        Class for RNA structure

        Input:
        -----------
        pdbfile: string
            PDB file with the RNA structure
        init_index: integer
            If provided the nulceotide are renumbered satrting from init_index.
            Default None (it keeps the nucleotide numbering of the PDB)
        test: boolean
            If True it runs a test for the class
            Default False
        """
        self.outdir = outdir
        self.test = test
        self.pdbfile = pdbfile
        self.pdbfile_name = pdbfile.split('/')[-1].split('.')[0]
        #print(self.pdbfile_name)
        self.biotite_pdb = pdb.PDBFile.read(self.pdbfile)
        self.biotite_atom_array = pdb.get_structure(self.biotite_pdb)[0]
        self.biotite_nucleotides = self.biotite_atom_array[
            struc.filter_nucleotides(self.biotite_atom_array)]
        self._get_nucleotides()
        # if init_index:
        #    self._reset_index(init_index)
        self._get_basepairs()
        self._get_glycosidic_bonds()
        self._get_edges()
        self._get_interactions()
        self.get_df_simple()

    # _reset_index NOT working!
    def _reset_index(self, init_index):
        """
        This function resets the index of the nucleotides to start from
        init_index, insted of the original index in the pdb file.

        Input variables:
        -----------
        init_index: integer
            Index used to reset the first nucleotide index.
        """
        new_index = init_index
        for i in range(len(self.biotite_nucleotides)):
            self.biotite_nucleotides[i] = new_index
            new_index += 1
        return 0

    def _get_nucleotides(self):
        """
        Computes a list of nucleotide ids and a list of nucleotide types which
        are stored as attributes, nucleotides_id and nucleotides_names,
        respectively.

        Returns:
            nucleotides_id: list (length number of nucleotides)
            nucleotides_names: list (length number of nucleotides)

        >>> rnaobj = RNA(pdbfile='test/trRosetta_1.pdb', test=True)
        >>> rnaobj._get_nucleotides()
        ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44], ['C', 'U', 'G', 'U', 'G', 'U', 'G', 'G', 'C', 'U', 'G', 'U', 'C', 'A', 'C', 'U', 'C', 'G', 'G', 'C', 'U', 'G', 'C', 'A', 'U', 'G', 'C', 'U', 'U', 'A', 'G', 'U', 'G', 'C', 'A', 'C', 'U', 'C', 'A', 'C', 'G', 'C', 'A', 'G'])
        """
        self.nucleotides_id = []
        self.nucleotides_names = []
        for nucleotide in struc.residue_iter(self.biotite_nucleotides):
            map_nucleotide, exact_match = struc.map_nucleotide(nucleotide)
            if map_nucleotide is None:
                continue
            self.nucleotides_id.append(nucleotide[0].res_id)
            if exact_match:
                self.nucleotides_names.append(map_nucleotide)
            else:
                self.nucleotides_names.append(map_nucleotide.lower())
        self.length = len(self.nucleotides_id)
        if self.test:
            return self.nucleotides_id, self.nucleotides_names
        else:
            return 0

    def _get_basepairs(self):
        """
        Computes a np.array of basepairs containing the index of the
        interacting nucleotides per basepair which is stored as an attribute
        (basepairs).

        Returns:
        -----------
            basepairs: list of length N where N is the number of
                       basepairs and where each element is a tupple with the
                       indexes of the two nucleotides forming the basepair

        >>> rnaobj = RNA(pdbfile='test/trRosetta_1.pdb', test=True)
        >>> rnaobj._get_basepairs()
        17
        """
        self.basepairs_atoms = struc.base_pairs(self.biotite_nucleotides)
        self.basepairs = struc.get_residue_positions(
            self.biotite_nucleotides, self.basepairs_atoms.flatten()).reshape(
            self.basepairs_atoms.shape)
        self.basepairs = [tuple(basepair) for basepair in self.basepairs]
        self.basepairs_id = [(self.nucleotides_id[basepair[0]],
                             self.nucleotides_id[basepair[1]])
                            for basepair in self.basepairs]
        if self.test:
            return len(self.basepairs)
        else:
            return 0

    def _get_glycosidic_bonds(self):
        """
        Computes the glycosidic interactions between the nitrogenous base and
        the sugar (cis or trans) for each nucleotide.

        Returns:
        -----------
        glycosidic_bonds: np array, shape=(N,)

            If the value of the orientation is:
            - 0: the glycosidic bond is in unknown conformation (INVALID).
            - 1: the glycosidic bond is in trans conformation.
            - 2: the glycosidic bond is in cis conformation.

        >>> rnaobj = RNA(pdbfile='test/trRosetta_1.pdb', test=True)
        >>> rnaobj._get_glycosidic_bonds()
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2]
        """
        self.glycosidic_bonds = struc.base_pairs_glycosidic_bond(
            self.biotite_nucleotides, self.basepairs_atoms)
        if self.test:
            return list(self.glycosidic_bonds)
        else:
            return 0

    def _get_edges(self):
        """
        Computes the interacting edge of each nucleotide in a basepair.

        Returns:
        -----------
        edges: np array, shape=(N, 2)

            If the value of the edge is:
            - 0: the nucleotide is not canonical or not HBonds present
                (INVALID).
            - 1: the nucleotide interacts with Watson-Crick edge.
            - 2: the nucleotide interacts with Hoogsteen edge.
            - 3: the nucleotide interacts with Sugar edge.

        Errors:
        -----------
        BadStructureError:
            If any edge is detected, the function raises an error.
            To avoid this, the function will return -1.


        >>> rnaobj = RNA(pdbfile='test/trRosetta_1.pdb', test=True)
        >>> rnaobj._get_edges()
        (17, 2)
        """
        try:
            self.edges = struc.base_pairs_edge(self.biotite_nucleotides,
                                               self.basepairs_atoms)
        except BadStructureError:
            self.edges = np.zeros((len(self.basepairs), 2), dtype=int)
            return -1
        if self.test:
            return self.edges.shape
        else:
            return 0

    def _get_interactions(self):
        """
        Assing to each nucleotide an specific type of RNA interaction

        Returns:
        -----------
        interactions: list of length number of nucleotides

            Example:
            ['cW', 'tW', 'cS', 'tS']

        >>> rnaobj = RNA(pdbfile='test/trRosetta_1.pdb', test=True)
        >>> rnaobj._get_interactions()
        ['cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'tW', 'tW']
        """

        interactions = []
        for bases, edges, orientation in zip(self.basepairs,
                                             self.edges,
                                             self.glycosidic_bonds):
            for base, edge in zip(bases, edges):
                if orientation == 1:
                    interaction = "c"
                elif orientation == 2:
                    interaction = "t"
                elif orientation == 0:
                    interaction = "x"
                if edge == 1:
                    interaction += "W"
                elif edge == 2:
                    interaction += "H"
                elif edge == 3:
                    interaction += "S"
                elif edge == 0:
                    interaction = "X"
                interactions.append(interaction)
        self.interactions = interactions
        if self.test:
            return self.interactions
        else:
            return 0

    def _get_hbonds(self):
        """
        Computes all hydrogen bonds in the RNA structure, not only the
        hydrogen bonds in the base pairs.

        Returns:
        -----------
        hbonds: numpy array, shape=(N, 3)
            Each row contains the indices of the atoms involved in
            the H-bond, with each column representing the Donor,
            Acceptor and Hydrogen.

            Example:
            Donor, Hydrogen, Acceptors
            [[566, 577,  14],
             [534, 546,  50],
             [505, 517,  79]]

        >>> rnaobj = RNA(pdbfile='test/trRosetta_1.pdb', test=True)
        >>> rnaobj._get_hbonds()
        (55, 3)
        """
        self.hbonds = struc.hbond(self.biotite_nucleotides)
        if self.test:
            return self.hbonds.shape
        else:
            return 0

    def get_df_simple(self):
        """
        Create a dataframe with the base pairs and the interaction type,
        without the hydrogen bonds.

        Returns:
        -----------
        df_simple: pandas dataframe
            Contains the base pairs and the interaction type.
            Example:
            BaseId1  BaseName1  BaseId2  BaseName2  Interaction Type
            1        G          18       C          cWcW
            2        A          17       U          tHcS

        >>> rnaobj = RNA(pdbfile='test/trRosetta_1.pdb', test=True)
        >>> rnaobj.get_df_simple()
        (17, 5)
        """
        df_simple = []

        for i in range(len(self.basepairs)):

            base1_index = self.basepairs[i][0]
            base2_index = self.basepairs[i][1]
            base1_id = self.nucleotides_id[base1_index]
            base2_id = self.nucleotides_id[base2_index]
            base1_name = self.nucleotides_names[base1_index]
            base2_name = self.nucleotides_names[base2_index]
            interaction_type1 = self.interactions[i*2]
            interaction_type2 = self.interactions[i*2+1]

            df_simple.append({
                'BaseId1': base1_id,
                'BaseName1': base1_name,
                'BaseId2': base2_id,
                'BaseName2': base2_name,
                'Interaction Type': interaction_type1 + interaction_type2
                })

        df_simple = pd.DataFrame(df_simple)
        self.df_simple = df_simple

        if self.test:
            return self.df_simple.shape
        else:
            return 0

    def _get_basepair_for_hbond(self, hbond):
        """
        Compute the basepair for a given hydrogen bond. The basepair is
        computed based on the index of the nucleotides involved in the
        hydrogen bond.

        Returns:
        -----------
        basepair: tuple
            Contains the basepair of the hydrogen bond with the donor
            and acceptor atoms. The hydrogen atom is not considered.
        """
        donor_index, hydrogen_index, acceptor_index = hbond
        donor_nucleotide = self.biotite_nucleotides[donor_index].res_id
        acceptor_nucleotide = self.biotite_nucleotides[acceptor_index].res_id
        basepair = (donor_nucleotide, acceptor_nucleotide)

        return basepair

    def _get_hbond_atom_names(self, hbond):
        """
        Compute the atom names for the donor and acceptor atoms in the
        hydrogen bond.

        """
        donor_name = self.biotite_nucleotides[hbond[0]].atom_name
        acceptor_name = self.biotite_nucleotides[hbond[2]].atom_name

        return donor_name, acceptor_name

    def _check_wobble(self, base1_name, base2_name, hbonds):
        """
        Given a basepair and a hydrogen bond, check if the interaction
        type is Wobble GU by checking the atoms involved in the hydrogen
        bond.

        Returns:
        -----------
        wobbleGU: integer
            Number of hydrogen bonds that are Wobble GU. If the number
            is 2 in the main code, the interaction type is Wobble GU.

        """

        wobbleGU = 0
        for hbond in hbonds:
            if hbond is not None:
                hbond = hbond.split('-')
                hbond_base1 = hbond[0]
                hbond_base2 = hbond[1]
                if base1_name == 'G' and base2_name == 'U':
                    if hbond_base1 == ('O6') and hbond_base2 == ('N3'):
                        wobbleGU += 1
                    elif hbond_base1 == ('N1') and hbond_base2 == ('O2'):
                        wobbleGU += 1

                elif base1_name == 'U' and base2_name == 'G':
                    if hbond_base1 == ('N3') and hbond_base2 == ('O6'):
                        wobbleGU += 1
                    elif hbond_base1 == ('O2') and hbond_base2 == ('N1'):
                        wobbleGU += 1

        if wobbleGU == 2:
            wobbleGU = True
        else:
            wobbleGU = False

        return wobbleGU

    def get_full_df(self):
        """
        Create a dataframe with the base pairs, the interaction type
        and the hydrogen bonds, including donors and acceptors atoms.
        Not only Watson-Crick, Hoogsteen and Sugar edges are considered,
        but also Wobble GU.

        Returns:
        -----------
        full_df: pandas dataframe
            Contains the base pairs, the interaction type and the
            hydrogen bonds.
            Example:

        BaseId1  BaseName1  BaseId2  BaseName2  InteractionType  HBond1  HBond2
        1        G          18       C          cWcW             O6-N4   N2-O2
        2        A          17       U          tHcS             N1-N3   N6-O4
        3        U          16       A          tHcW             N3-O2

        >>> rnaobj = RNA(pdbfile='test/trRosetta_1.pdb', test=True)
        >>> rnaobj.get_full_df()
        ((17, 8), 'WobbleGU')
        """
        if not hasattr(self, 'hbonds'):
            self._get_hbonds()

        df_complex = self.df_simple.copy()

        # Compute a dictionary which given a basepair return its hbonds
        basepair_hbonds = {}
        basepair_hbonds_names = {}
        for hbond in self.hbonds:
            basepair = self._get_basepair_for_hbond(hbond)
            donor_name, acceptor_name = self._get_hbond_atom_names(hbond)
            if basepair not in self.basepairs_id:
                _basepair = (basepair[1], basepair[0])
                if _basepair in self.basepairs_id:
                    basepair = (_basepair[0], _basepair[1])
                    hbond_name = '%s-%s' % (acceptor_name, donor_name)
                else:
                    continue
            else:
                basepair = (basepair[0], basepair[1])
                hbond_name = '%s-%s' % (donor_name, acceptor_name)
            if basepair not in basepair_hbonds.keys():
                basepair_hbonds[basepair] = [hbond]
                basepair_hbonds_names[basepair] = [hbond_name]
            else:
                basepair_hbonds[basepair].append(hbond)
                basepair_hbonds_names[basepair].append(hbond_name)

        # Fill df_complex with the new HBonds
        self.max_hbonds = max([len(hbonds) for hbonds in
                        basepair_hbonds_names.values()])
        for i in range(self.max_hbonds):
            df_complex['Hbond' + str(i+1)] = None
        for basepair, hbond_names in basepair_hbonds_names.items():
            for i, hbond_name in enumerate(hbond_names):
                df_complex.loc[(df_complex['BaseId1'] == basepair[0])
                             & (df_complex['BaseId2'] == basepair[1]),
                             'Hbond' + str(i+1)] = hbond_name

        # Check if interaction type is Wobble GU using hbonds information
        for index, row in df_complex.iterrows():
            hbonds = []
            base1_name = row['BaseName1']
            base2_name = row['BaseName2']
            for i in range(self.max_hbonds):
                hbond = row['Hbond' + str(i+1)]
                hbonds.append(hbond)
            wobbleGU = self._check_wobble(base1_name, base2_name, hbonds)
            if wobbleGU:
                df_complex.loc[index, 'Interaction Type'] = 'WobbleGU'

        self.df_complex = df_complex
        #print(self.df_complex)
        if self.test:
            #return self.df_complex.shape
            return self.df_complex.shape, self.df_complex.loc[14, 'Interaction Type']
        else:
            return 0

    def save_df_simple(self, outname):
        """
        Store the dataframe in a txt/csv file. The df stored
        will be the one the user has computed, simple or complex.
        """

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        out_file = os.path.join(outdir, f'dataframe.txt')
        df.to_csv(out_file, type='txt', index=False, sep='\t')

    def save_df_complex(self, outname):
        """

        """
        # has de mirar si df_complex existeix com attribute  del teu objte
        # en cas contrari calcular-lo abans de guardar
        return

    def plot_df_interactions(self, merged_df):
        """
        Plot the number of interactions per base pair with a list of all
        base pairs in the x-axis and the number of interactions in the y-axis.
        Given multiple files, the user will be able to compare the interactions
        between different structures, as each interaction will be represented
        with a different color.

        Returns:
        -----------
        interactions.png: png file
            Plot with the number of interactions per base pair.

        # >>> rnaobj = RNA(pdbfile='test/trRosetta_1.pdb', test=True, outdir='test')
        # >>> rnaobj.plot_df_interactions(merged_df=pd.DataFrame())
        # interactions.png
        """
        # Test is not working for plot_df_interactions
        pair_df = []
        pair_df_index = {}
        self.merged_df = merged_df
        for index, row in self.merged_df.iterrows():
            bp_label = f"{row['BaseName1']}{row['BaseId1']}-{row['BaseName2']}{row['BaseId2']}"
            if bp_label not in pair_df_index:
                pair_df_index[bp_label] = len(pair_df)
                pair_df.append(bp_label)

        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(28, 15))
        x = np.arange(len(pair_df))
        width = 0.35
        unique_interactions = self.merged_df['Interaction Type'].unique()
        color_palette = sns.color_palette("tab10", len(unique_interactions))

        total_interactions_per_pair = np.zeros((len(pair_df), len(unique_interactions)))
        for i, pair in enumerate(pair_df):
            total_interactions = []
            for j, interaction in enumerate(unique_interactions):
                count = self.merged_df[(self.merged_df['BaseName1'] + self.merged_df['BaseId1'].astype(str) + "-" + self.merged_df['BaseName2'] + self.merged_df['BaseId2'].astype(str) == pair) &
                                   (self.merged_df['Interaction Type'] == interaction)].shape[0]
                total_interactions_per_pair[i, j] = count

        bottom = None
        for i, interaction in enumerate(unique_interactions):
            plt.bar(x, total_interactions_per_pair[:, i], width, label=interaction, bottom=bottom, color=color_palette[i])
            if bottom is None:
                bottom = total_interactions_per_pair[:, i]
            else:
                bottom += total_interactions_per_pair[:, i]

        self.plot_name = f'interactiom_barplot.png'
        plt.xlabel('Base Pair', fontsize=22)
        plt.ylabel('Count', fontsize=22)
        plt.title('Number of interactions per base pair', fontsize=26)
        plt.xticks(x, pair_df, fontsize=14, rotation=45)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=16)
        if self.outdir:
            plt.savefig(os.path.join(self.outdir, self.plot_name))
        else:
            plt.savefig(self.plot_name)
        plt.show()

        if self.test:
            return self.plot_name
        else:
            return 0

    def merge_df(self, merged_df):
        """
        Merge the dataframes of the different structures to generate a
        unique dataframe with all the interactions.

        Input:
        -----------
        df: pandas dataframe
            Dataframe to be merged.

        """

        rnaobj.get_full_df()
        merged_df = pd.concat([merged_df, self.df_complex])

        return merged_df

"""
--------------------------------------
--------------------------------------
--------------------------------------
"""

def get_files(directory):
    """
    Get the files in the directory.

    Input:
    -----------
    directory: string
        Path to the directory containing the pdb files.

    Returns:
    -----------
    files: list of strings
        List of the pdb files in the directory.

    """
    files = glob.glob(os.path.join(directory, '*.pdb'))
    return files


if __name__ == '__main__':
    import doctest
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--pdb',
                        help='Provide the individual file names to the structures '+
                        'to be analyzed',)
    parser.add_argument('--test',
                        help='Test the code',
                        action='store_true')
    parser.add_argument('--dir',
                        help='Provide the directory with multiple pdb files to be analyzed')
    parser.add_argument('--dataframe',
                        help='Choose the dataframe to be stored: Simple or Complex. Default=simple',
                        choices=['simple', 'complex'],
                        default='simple')
    parser.add_argument('--out',
                        help='Path to the output directory')

    args = parser.parse_args()

    if args.test:
        doctest.testmod(
            optionflags=doctest.ELLIPSIS | doctest.REPORT_ONLY_FIRST_FAILURE)
        sys.exit()

    if args.pdb:
        #print(args.pdb)
        rnaobj = RNA(pdbfile=args.pdb)
        if args.dataframe == 'complex':
            merged_df = pd.DataFrame()
            merged_df = rnaobj.merge_df(merged_df)
            rnaobj.plot_df_interactions(merged_df)

    if args.dir:
        files = get_files(args.dir)
        merged_df = pd.DataFrame()
        for file in files:
            #print(file)
            rnaobj = RNA(pdbfile=file)
            if args.dataframe == 'complex':
                merged_df = rnaobj.merge_df(merged_df)
        rnaobj.plot_df_interactions(merged_df)


    # print(rnaobj.df_interactions)

