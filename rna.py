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
    def __init__(self, pdbfile, init_index=None, test=False):
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
        self.test = test
        self.biotite_pdb = pdb.PDBFile.read(pdbfile)
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

        >>> rnaobj = RNA(pdbfile='test/trRosetta_2.pdb', test=True)
        >>> rnaobj._get_nucleotides()
        ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], ['G', 'A', 'U', 'C', 'U', 'C', 'U', 'U', 'G', 'U', 'A', 'G', 'A', 'U', 'C'])
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

        Return:
            basepairs: list of length N where N is the number of
                       basepairs and where each element is a tupple with the
                       indexes of the two nucleotides forming the basepair

        >>> rnaobj = RNA(pdbfile='test/trRosetta_2.pdb', test=True)
        >>> rnaobj._get_basepairs()
        5
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

        >>> rnaobj = RNA(pdbfile='test/trRosetta_2.pdb', test=True)
        >>> rnaobj._get_glycosidic_bonds()
        [1, 1, 1, 1, 1]
        """
        self.glycosidic_bonds = struc.base_pairs_glycosidic_bond(
            self.biotite_nucleotides, self.basepairs_atoms)
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


        >>> rnaobj = RNA(pdbfile='test/trRosetta_2.pdb', test=True)
        >>> rnaobj._get_edges()
        (5, 2)
        """
        try:
            self.edges = struc.base_pairs_edge(
                self.biotite_nucleotides, self.basepairs_atoms)
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

        >>> rnaobj = RNA(pdbfile='test/trRosetta_2.pdb', test=True)
        >>> rnaobj._get_interactions()
        ['cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW', 'cW']
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
            [[566, 577,  14],
             [534, 546,  50],
             [505, 517,  79]]

        >>> rnaobj = RNA(pdbfile='test/trRosetta_2.pdb', test=True)
        >>> rnaobj._get_hbonds()
        (22, 3)
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

        >>> rnaobj = RNA(pdbfile='test/trRosetta_2.pdb', test=True)
        >>> rnaobj.get_df_simple()
        (5, 5)
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

        """
        donor_index, hydrogen_index, acceptor_index = hbond
        donor_nucleotide = self.biotite_nucleotides[donor_index].res_id
        acceptor_nucleotide = self.biotite_nucleotides[acceptor_index].res_id
        basepair = (donor_nucleotide, acceptor_nucleotide)

        return basepair

    def _get_hbond_atom_names(self, hbond):
        """

        """
        donor_name = self.biotite_nucleotides[hbond[0]].atom_name
        acceptor_name = self.biotite_nucleotides[hbond[2]].atom_name

        return donor_name, acceptor_name

    def _check_wobble(self, base1_name, base2_name, hbond):
        """

        """
        hbond_base1 = hbond[0]
        hbond_base2 = hbond[1]
        wobbleGU = False
        if base1_name == 'G' and base2_name == 'U':
            if hbond_base1 == ('O6' or 'N1') and hbond_base2 == ('N3' or 'O4'):
                wobbleGU = True

        elif base1_name == 'U' and base2_name == 'G':
            if hbond_base1 == ('O4' or 'N3') and hbond_base2 == ('O6' or 'N1'):
                wobbleGU = True

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
            hbond_name = '%s-%s' % (donor_name, acceptor_name)
            if basepair not in self.basepairs_id:
                _basepair = (basepair[1], basepair[0])
                if _basepair in self.basepairs_id:
                    basepair = (_basepair[0], _basepair[1])
                else:
                    continue
            else:
                basepair = (basepair[0], basepair[1])
            if basepair not in basepair_hbonds.keys():
                basepair_hbonds[basepair] = [hbond]
                basepair_hbonds_names[basepair] = [hbond_name]
            else:
                basepair_hbonds[basepair].append(hbond)
                basepair_hbonds_names[basepair].append(hbond_name)

        # Fill df_complex with the new HBonds
        max_hbonds = max([len(hbonds) for hbonds in
                        basepair_hbonds_names.values()])
        for i in range(max_hbonds):
            df_complex['Hbond' + str(i+1)] = None

        for basepair, hbond_names in basepair_hbonds_names.items():
            for i, hbond_name in enumerate(hbond_names):
                df_complex.loc[(df_complex['BaseId1'] == basepair[0])
                             & (df_complex['BaseId2'] == basepair[1]),
                             'Hbond' + str(i+1)] = hbond_name

        # Check if interaction type is Wobble GU using hbonds information
        for index, row in self.df_simple.iterrows():
            basepair = (row['BaseId1'], row['BaseId2'])
            hbonds = basepair_hbonds[basepair]
            for hbond in hbonds:
                wobbleGU = self._check_wobble(row['BaseName1'],
                                              row['BaseName2'],
                                              hbond)
                if wobbleGU:
                    df_complex['InteractionType'] = 'Wobble GU'
        self.df_complex = df_complex

        return 0

    def store_df(self):
        """
        """

        global_df = pd.concat([set.full_df for set in self.sets])

        self.global_df = global_df

        return 0

    def plot_df_interactions(self, indir, outdir):
        """
        """
        files = glob.glob(os.path.join(indir, '*.pdb'))
        global_df = pd.DataFrame

        for file in files:
            rnaobj = RNA(pdbfile=file)
            rnaobj.get_main() # Check whether it is correctly implemented
            global_df = pd.concat([global_df, rnaobj.full_df])

        pair_df = []
        pair_df_index = {}

        for index, row in global_df.iterrows():
            bp_label = f"{row['BaseName1']}{row['BaseId1']}-{row['BaseName2']}{row['BaseId2']}"
            if bp_label not in pair_df_index:
                pair_df_index[bp_label] = len(pair_df)
                pair_df.append(bp_label)

        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(28, 15))
        x = np.arange(len(pair_df))
        width = 0.35
        unique_interactions = global_df['InteractionType'].unique()
        color_palette = sns.color_palette("tab10", len(unique_interactions))

        total_interactions_per_pair = np.zeros((len(pair_df), len(unique_interactions)))
        for i, pair in enumerate(pair_df):
            total_interactions = []
            for j, interaction in enumerate(unique_interactions):
                count = global_df[(global_df['BaseName1'] + global_df['BaseId1'].astype(str) + "-" + global_df['BaseName2'] + global_df['BaseId2'].astype(str) == pair) &
                                   (global_df['InteractionType'] == interaction)].shape[0]
                total_interactions_per_pair[i, j] = count

        bottom = None
        for i, interactin in enumerate(unique_interactions):
            plt.bar(x, total_interactions_per_pair[:, i], width, label=interaction, bottom=bottom, color=color_palette[i])
            if bottom is None:
                bottom = total_interactions_per_pair[:, i]
            else:
                bottom += total_interactions_per_pair[:, i]

        plt.xlabel('Base Pair', fontsize=22)
        plt.ylabel('Count', fontsize=22)
        plt.title('Number of interactions per base pair', fontsize=26)
        plt.xticks(x, pair_df, fontsize=14, rotation=45)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=16)
        plt.savefig(os.path.join(outdir, 'interactions.pdf'))
        plt.show()

        return 0

# Args to be implemented:
# --dir: path to the directory containing the pdb files
# --output: Simple or Complex

if __name__ == '__main__':
    import doctest
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--pdb',
                        help='',)
    parser.add_argument('--test',
                        help='Test the code',
                        action='store_true')
    args = parser.parse_args()

    if args.test:
        doctest.testmod(
            optionflags=doctest.ELLIPSIS | doctest.REPORT_ONLY_FIRST_FAILURE)
        sys.exit()

    rnaobj = RNA(pdbfile=args.pdb)

    rnaobj.get_full_df()
    # print(rnaobj.df_interactions)

    # rnaobj.plot_df_interactions(outname='.pdf')
