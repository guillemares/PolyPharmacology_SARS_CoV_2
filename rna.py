#!/usr/bin/env python
import biotite.structure.io.pdb as pdb
import biotite.structure as struc
from biotite.structure.error import BadStructureError
import sys
import numpy as np


class RNA(object):
    """
    >>> RNA(pdbfile='test/3dRNA_109.min.pdb')
    """
    def __init__(self, pdbfile, init_index=None):
        """

        """
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
        self.get_df_interactions()

    def _reset_index(self, init_index):
        """

        """
        # self.nucleotides_id =
        return 0

    def _get_nucleotides(self):
        """
        Parameters:
        -----------
        biotite_atom_array: AtomArray
            Contains the structure of a nucleic acid.

        Returns:
        -----------
        nucleotides_id: list
            Contains the nucleotides numbers (id) in the structure.
        nucleotides_names: list
            Contains the nucleotides names in the structure.
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

    def _get_basepairs(self):
        """
        Parameters:
        -----------
        biotite_atom_array: AtomArray
            Contains the structure of a nucleic acid.

        Returns:
        -----------
        basepairs: numpy array, shape=(N, 2)
            Contains the base pairs in the structure. Each row contains the
            residue indices of a base pair.

            Example:
            [[1, 8], [2, 7], [3, 6], [4, 5]]

        """
        self.basepairs_atoms = struc.base_pairs(self.biotite_nucleotides)   # Gives atoms index. It just returns one index per base pair, so it is not useful
        self.basepairs = struc.get_residue_positions(                       # for hydrogen bonding analysis, only for base pairing.
            self.biotite_nucleotides, self.basepairs_atoms.flatten()
        ).reshape(self.basepairs_atoms.shape)                               # Gives residue index
        return 0

    def _get_glycosidic_bonds(self):
        """
        Parameters:
        -----------
        biotite_atom_array: AtomArray
            Contains the structure of a nucleic acid.
        basepairs: numpy array, shape=(N, 2)
            Contains the base pairs in the structure.

        Returns:
        -----------
        glycosidic_bonds: numpy array, shape=(N,)
            Contains the glycosidic bonds in the structure. Each element.
            contains the glycosidic bond of a base pair with the shape (N,).

            If the value of the orientation is:
            - 0: the glycosidic bond is in unknown conformation (INVALID).
            - 1: the glycosidic bond is in trans conformation.
            - 2: the glycosidic bond is in cis conformation.

            Example:
            [1, 2, 1, 2]

        """
        self.glycosidic_bonds = struc.base_pairs_glycosidic_bond(
            self.biotite_nucleotides, self.basepairs)
        return 0

    def _get_edges(self):
        """
        Parameters:
        -----------
        atom_array: AtomArray
            Contains the structure of a nucleic acid.
        basepairs: numpy array, shape=(N, 2)
            Contains the base pairs in the structure.

        Returns:
        -----------
        edges: numpy array, shape=(N, 2)
            Contains the edges in the structure. 

            If the value of the edge is:
            - 0: the nucleotide is not canonical or not HBonds present (INVALID).
            - 1: the nucleotide interacts with Watson-Crick edge.
            - 2: the nucleotide interacts with Hoogsteen edge.
            - 3: the nucleotide interacts with Sugar edge.

            Example:
            [[1, 1], [2, 2], [3, 2], [1, 3]]

        Errors:
        -----------
        BadStructureError:
            If any edge is detected, the function raises an error.
            To avoid this, the function will return -1.
        """
        try:
            self.edges = struc.base_pairs_edge(
                self.biotite_nucleotides, self.basepairs)
        except BadStructureError:
            # self.edges = None
            self.edges = np.zeros(self.basepairs.shape, dtype=int)
            return -1
        return 0

    def _get_interactions(self):
        """
        Parameters:
        -----------
        basepairs: numpy array, shape=(N, 2)
            Contains the base pairs in the structure.
        edges: numpy array, shape=(N, 2)
            Contains the edges in the structure.
        glycosidic_bonds: numpy array, shape=(N,)
            Contains the glycosidic bonds in the structure.

        Returns:
        -----------
        interactions: list of str
            Contains the interactions between the base pairs in the structure.
            Example:
            ['cW', 'tW', 'cS', 'tS']

        """

        interactions = []
        for bases, edges, orientation in zip(
            self.basepairs, self.edges, self.glycosidic_bonds):
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
        return 0
    
    def _get_base_pairs_interactions(self):
        """
        
        Adjusts the interactions between the base pairs in the structure with the
        correct nucleotides names (and correct the index to start from 1 instead of 0).

        Parameters:
        -----------
        basepairs: numpy array, shape=(N, 2)
            Contains the base pairs in the structure.
        interactions: list
            Contains the interactions between the base pairs in the structure.
        nucleotides_names: list
            Contains the nucleotides names in the structure.

        Returns:
        -----------
        base_pairs_interactions: list of list
            Contains the interactions between the base pairs in the structure
            with the nucleotides names.
            Each element contains the base pair index, the nucleotide name and the
            interaction. 
            Example:
            [[[1, 'A', 'cW'], [2, 'U', 'tW']], [[3, 'G', 'cS'], [4, 'C', 'tS']]]

        """
        base_pairs_interactions = []
        for i in range(self.basepairs.shape[0]):
            edge1 = self.interactions[i*2]
            edge2 = self.interactions[i*2+1]

            base_pairs1 = self.basepairs[i, 0] + 1
            base_pairs2 = self.basepairs[i, 1] + 1

            base_pairs_interactions.append(
                [base_pairs1, self.nucleotides_names[self.basepairs[i, 0]], edge1],
                [base_pairs2, self.nucleotides_names[self.basepairs[i, 1]], edge2]
            )

        self.base_pairs_interactions = base_pairs_interactions
        return 0

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

    def _get_hbonds_for_base_pair(self):
        """
        Parameters:
        """
        self.hbonds_for_base_pair = {}
        for hbond in self.hbonds:
            donor_index = hbond[0]
            hydrogen_index = hbond[1]
            acceptor_index = hbond[2]

            donor_nucleotide = self.biotite_nucleotides[donor_index].res_id
            hydrogen_nucleotide = self.biotite_nucleotides[hydrogen_index].res_id
            acceptor_nucleotide = self.biotite_nucleotides[acceptor_index].res_id

            for base_pair in self.basepairs:
                if (donor_nucleotide in base_pair[0] and acceptor_nucleotide in base_pair[1]) \
                or (donor_nucleotide in base_pair[1] and acceptor_nucleotide in base_pair[0]):
                    base_pair_key = (tuple(base_pair[0]), tuple(base_pair[1]))
                    if base_pair_key not in self.hbonds_for_base_pair:
                        self.hbonds_for_base_pair[base_pair_key] = []
                    
                    donor_atom_name = self.nucleotides_names[donor_nucleotide].atom_name       # Returns the donor atom type
                    hydrogen_atom_name = self.nucleotides_names[hydrogen_nucleotide].atom_name # Returns the hydrogen atom type
                    acceptor_atom_name = self.nucleotides_names[acceptor_nucleotide].atom_name # Returns the acceptor atom type

                    if donor_nucleotide == base_pair[0][0]:
                        base1_atom = donor_atom_name


        return 0

    def get_df_interactions(self):
        """

        """
        self._get_interactions()
        self._get_hbonds()
        self._get_hbonds_for_base_pair()
        # self.df_interactions
        return 0

    def plot_df_interactions(self, outname):
        return 0


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
