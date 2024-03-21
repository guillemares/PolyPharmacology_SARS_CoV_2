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
        if init_index:
            self._reset_index(init_index)
        self._get_basepairs()
        self._get_glycosidic_bonds()
        self._get_edges()
        self._get_interactions()
        # self._get_hbonds_for_basepair()
        self.get_simple_df() # This function will be used only when the user types --output = 'Simple'
        # self.get_full_df() # This function will be used only when the user types --output = 'Complex'
        # self.plot_df_interactions() out of the loop, it has to be used after storing all DataFrames in global_df

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
            basepairs: np array, shape=(N,2) where N is the number of
                       basepairs

        >>> rnaobj = RNA(pdbfile='test/trRosetta_2.pdb', test=True)
        >>> rnaobj._get_basepairs()
        (5, 2)
        """
        self.basepairs_atoms = struc.base_pairs(self.biotite_nucleotides)   # Gives atoms index. It just returns one index per base pair, so it is not useful
        self.basepairs = struc.get_residue_positions(
            self.biotite_nucleotides, self.basepairs_atoms.flatten()).reshape(
            self.basepairs_atoms.shape)
        if self.test:
            return self.basepairs.shape
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
            self.edges = np.zeros(self.basepairs.shape, dtype=int)
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

        # def _get_basepairs_interactions(self):
        #     """
        #     Assign to each basepair a 2 dimensional np array where each
        #     element of it containts the following information:
        #     ['nucleotide1 id', 'nucleotide1 name', 'interaction type']

        #     Returns:
        #     -----------
        #     basepairs_interactions: np array shape (N,2,3)

        #     Example:
        #     [[[1, 'A', 'cW'], [2, 'U', 'tW']], [[3, 'G', 'cS'], [4, 'C', 'tS']]]

        #     >>> rnaobj = RNA(pdbfile='test/trRosetta_2.pdb', test=True)
        #     >>> rnaobj._get_basepairs_interactions()
        #     (5, 2, 3)
        #     """
        #     basepairs_interactions = []
        #     for i in range(self.basepairs.shape[0]):
        #         edge1 = self.interactions[i*2]
        #         edge2 = self.interactions[i*2+1]

        #         basepairs1 = self.basepairs[i, 0] + 1
        #         basepairs2 = self.basepairs[i, 1] + 1

        #         basepairs_interactions.append(
        #             [[basepairs1, self.nucleotides_names[self.basepairs[i, 0]],
        #              edge1],
        #              [basepairs2, self.nucleotides_names[self.basepairs[i, 1]],
        #              edge2]])

        #     self.basepairs_interactions = np.asarray(basepairs_interactions)
        #     if self.test:
        #         return self.basepairs_interactions.shape
        #     else:
        #         return 0

    def _get_hbonds(self):
        """
        Parameters:
        -----------
        atom_array: AtomArray
            Contains the structure of a nucleic acid.

        Returns:
        -----------
        hbonds: numpy array, shape=(N, 3)
            Contains the hydrogen bonds in the structure. Each row contains
            the indices of the atoms involved in a hydrogen bond, with the
            columns representing the Donor, Acceptor and Hydrogen atom.

            Example:
            [[566, 577,  14],
             [534, 546,  50],
             [505, 517,  79]]

        """
        self.hbonds = struc.hbond(self.biotite_nucleotides)
        return 0

    def _get_hbonds_for_basepair(self):
        """
        Associate the hydrogen bonds with the base pairs in the structure,
        using the index of the atoms (Donor, Acceptor and Hydrogen).
        Note that each base pair may have more than one hydrogen bond.

        Parameters:
        -----------
        hbonds: numpy array, shape=(N, 3)
            Contains the hydrogen bonds in the structure. Each row contains
            the indices of the atoms involved in a hydrogen bond, with the
            columns representing the Donor, Acceptor and Hydrogen atom.
        basepairs: numpy array, shape=(N, 2)
            Contains the base pairs in the structure.

        Returns:
        -----------
        hbonds_for_basepair: dict
            Contains the hydrogen bonds for each base pair in the structure.
            The key is a tuple with the base pair, and the value is a list
            with the atoms involved in the hydrogen bond (Donor, Acceptor).
            Hydrogen atom is not included for simplicity.

            Example:
            {((1, 'G', 'cW'), (18, 'C', 'cW')): [('O6', 'N4'), ('N2', 'O2'), ('N1', 'N3')],
             ((2, 'A', 'cW'), (17, 'U', 'cW')): [('N1', 'N3'), ('N6', 'O4')]}
        """
        self.hbonds_for_basepair = {}
        for hbond in self.hbonds:
            donor_index = hbond[0]
            hydrogen_index = hbond[1]
            acceptor_index = hbond[2]

            donor_nucleotide = self.biotite_nucleotides[donor_index].res_id
            hydrogen_nucleotide = self.biotite_nucleotides[hydrogen_index].res_id
            acceptor_nucleotide = self.biotite_nucleotides[acceptor_index].res_id

            for basepair in self.basepairs:
                if (donor_nucleotide in basepair[0] and acceptor_nucleotide in basepair[1]) \
                or (donor_nucleotide in basepair[1] and acceptor_nucleotide in basepair[0]):
                    basepair_key = (tuple(basepair[0]), tuple(basepair[1]))
                    if basepair_key not in self.hbonds_for_basepair:
                        self.hbonds_for_basepair[basepair_key] = []
                    
                    donor_atom_name = self.nucleotides_names[donor_nucleotide].atom_name        # Returns the donor atom type
                    hydrogen_atom_name = self.nucleotides_names[hydrogen_nucleotide].atom_name  # Returns the hydrogen atom type
                    acceptor_atom_name = self.nucleotides_names[acceptor_nucleotide].atom_name  # Returns the acceptor atom type

                    if donor_nucleotide == basepair[0][0]:
                        base1_atom = donor_atom_name
                        base2_atom = acceptor_atom_name
                    else:
                        base1_atom = acceptor_atom_name
                        base2_atom = donor_atom_name

                    self.hbonds_for_basepair[basepair_key].append([base1_atom, base2_atom])

        return 0

    # I propose two possibilities. It is similar to my original code but implementing the idea 
    # that the user may prefer a simpler calculation. Arguments --output = (Simple or Complex).
    # a) If the user just needs to compute de base pairs a simpler dataframe
    #    can be created with these information:
    #    'BaseId1', 'BaseName1', 'BaseId2', 'BaseName2', 'Annotation' (annotation 1 and 2)
    #
    # b) If the user needs to compute the hydrogen bonds, a more complex dataframe
    #    can be created with these information:
    #    'BaseId1', 'BaseName1', 'BaseId2', 'BaseName2', 'Annotation', 'HBond1', 'HBond2', 'HBond3', ..., 'WobbleGU?'

    def get_simple_df(self):
        """
        Create a dataframe with the base pairs and the interaction type, 
        without the hydrogen bonds. 
        
        Parameters:
        -----------
        basepairs: numpy array, shape=(N, 2)
            Contains the base pairs in the structure.
        nucleotides_names: list
            Contains the nucleotides names in the structure.

        Returns:
        -----------
        simple_df: pandas dataframe
            Contains the base pairs and the interaction type.
            Example:
            BaseId1  BaseName1  BaseId2  BaseName2  Interaction Type
            1        G          18       C          cWcW
            2        A          17       U          tHcS

        """
        simple_df = []

        for basepair in self.basepairs:
            base1_id = basepair[0][0]
            base1_name = basepair[0][1]
            interaction_type1 = basepair[0][2]

            base2_id = basepair[1][0]
            base2_name = basepair[1][1]
            interaction_type2 = basepair[1][2]
            
            simple_df.append({
                'BaseId1': base1_id,
                'BaseName1': base1_name,
                'BaseId2': base2_id,
                'BaseName2': base2_name,
                'Interaction Type': interaction_type1 + interaction_type2})

        simple_df = pd.DataFrame(simple_df)
        self.simple_df = simple_df

        return 0

    def get_full_df(self):
        """
        Create a dataframe with the base pairs, the interaction type and the hydrogen bonds,
        including donors and acceptors atoms. Not only Watson-Crick, Hoogsteen and Sugar edges
        are considered, but also Wobble GU.

        Parameters:
        -----------
        basepairs: numpy array, shape=(N, 2)
            Contains the base pairs in the structure.
        nucleotides_names: list
            Contains the nucleotides names in the structure.
        hbonds_for_basepair: dict
            Contains the hydrogen bonds for each base pair in the structure.
            The key is a tuple with the base pair, and the value is a list
            with the atoms involved in the hydrogen bond (Donor, Acceptor).
            Hydrogen atom is not included for simplicity.

        Returns:
        -----------
        full_df: pandas dataframe
            Contains the base pairs, the interaction type and the hydrogen bonds.
            Example:
            BaseId1  BaseName1  BaseId2  BaseName2  Interaction Type  HBond1  HBond2  HBond3
            1        G          18       C          cWcW              O6-N4   N2-O2   N1-N3
            2        A          17       U          tHcS              N1-N3   N6-O4        
            3        U          16       A          tHcW              N3-O2                
        """
        if not hasattr(self, 'hbonds'):
            self._get_hbonds()

        full_df = []

        for basepair, self.hbonds in self.hbonds_for_basepair.items():
            base1_id = basepair[0][0]
            base1_name = basepair[0][1]
            interaction_type1 = basepair[0][2]
            
            base2_id = basepair[1][0]
            base2_name = basepair[1][1]
            interaction_type2 = basepair[1][2]
            
            # no need
            data_full_df = {
                'BaseId1': base1_id,
                'BaseName1': base1_name,
                'BaseId2': base2_id,
                'BaseName2': base2_name
            }

            i = 1

            has_wobble = False

            for hbond in self.hbonds:
                hbond_base1 = hbond[0]
                hbond_base2 = hbond[1]
                wobbleGU = False

            # create function to detect if wobble or not _check_wobble(base1,base2,atom1,atom2)
                if base1_name == 'G' and base2_name == 'U':
                    if hbond_base1 ==  ('O6' or 'N1') and hbond_base2 == ('N3' or 'O4'):
                        wobbleGU = True
                
                elif base1_name == 'U' and base2_name == 'G':
                    if hbond_base1 ==  ('O4' or 'N3') and hbond_base2 == ('O6' or 'N1'):
                        wobbleGU = True
            
                has_wobble = has_wobble or wobbleGU
                data_full_df[f'HBond{i}'] = f'{hbond_base1}-{hbond_base2}'
                
                i += 1

            if has_wobble:
                data_full_df['InteractionType'] = 'Wobble GU'
            else:
                data_full_df['InteractionType'] = interaction_type1 + interaction_type2

            full_df.append(data_full_df)

        full_df = pd.DataFrame(full_df)
        column_order = ['BaseId1', 'BaseName1', 'BaseId2', 'BaseName2']
        hbond_columns = [column for column in full_df.columns if column.startswith('HBond')]
        interaction_column = ['InteractionType']
        new_columns = column_order + interaction_column + hbond_columns
        full_df = full_df[new_columns]
        self.full_df = full_df

        return 0


    def get_main(self):
        """

        """
        self._get_interactions()
        self._get_hbonds()
        self._get_hbonds_for_basepair()
        self._get_basepairs_interactions()
        self.get_simple_df()
        self.get_full_df()
        self.plot_df_interactions() # !!!!

        return 0

    def store_df(self):
        """
        """

        global_df = pd.concat([set.full_df for set in self.sets])
        
        self.global_df = global_df
        
        return 0

    # As it needs multiple files to be created, it should be
    # used when the user types --dir = 'path/to/dir', where
    # all pdb files will be read.

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

    # print(rnaobj.df_interactions)

    # rnaobj.plot_df_interactions(outname='.pdf')
