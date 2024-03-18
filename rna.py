#!/usr/bin/env python
import biotite.structure.io.pdb as pdb
import biotite.structure as struc
from biotite.structure.error import BadStructureError


class RNA(object):
    """

    """
    def __init__(self, pdbfile, init_index=None):
        """

        """
        self.biotite_pdb = pdb.PDBFile.read(pdbfile)
        self.biotite_atom_array = pdb.get_structure(self.biotite_pdb)[0]
        self.biotite_nucleotides = self.biotite_atom_array[struc.filter_nucleotides(self.biotite_atom_array)]
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
        quina info utilitza
        que torna com ho torna
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

        """
        self.basepairs = struc.base_pairs(self.biotite_nucleotides)
        # self.basepairs_flatten = 
        return 0

    def _get_glycosidic_bonds(self):
        """
        return a list of .... where 1 is cys and 2 is trans
        """
        self.glycosidic_bonds = struc.base_pairs_glycosidic_bond(self.biotite_nucleotides, self.basepairs)
        return 0

    def _get_edges(self):
        """
        return a numpy array .... wehre 0 means ... where 1 means .. 2 means ..
        """
        try:
            self.edges = struc.base_pairs_edge(self.biotite_nucleotides, self.basepairs)
        except BadStructureError:
            self.edges = None
            return -1
        return 0

    def _get_interactions(self):
        """
        self.basepairs_flatten
        """

        interactions = []
        self.interactions = interactions
        return 0

    def _get_hbonds(self):
        """

        """
        self.hbonds = struc.hbond(self.biotite_nucleotides)
        return 0

    def _get_hbonds_for_base_pair(self):
        """

        """
        self.hbonds_for_base_pair = {}

        return 0

    def get_df_interactions(self):
        """

        """
        self._get_interactions()
        self._get_hbonds()
        self._get_hbonds_for_base_pair()
        self.df_interactions
        return 0

    def plot_df_interactions(self, outname):
        return 0

if __name__ == '__main__':
    print('hello')
    pdbfile = '3dRNA_109.min.pdb'
    rnaobj = RNA(pdbfile=pdbfile)
    print(rnaobj.df_interactions)
    rnaobj.plot_df_interactions(outname='.pdf')
