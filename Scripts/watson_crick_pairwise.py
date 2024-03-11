from tempfile import gettempdir
import biotite
import biotite.structure.io.pdb as pdb
import biotite.database.rcsb as rcsb
import biotite.structure as struc
import numpy as np

pdb_file_path = rcsb.fetch("6ZYB", "pdb", gettempdir())
pdb_file = pdb.PDBFile.read(pdb_file_path)
atom_array = pdb.get_structure(pdb_file)[0]
nucleotides = atom_array[struc.filter_nucleotides(atom_array)]

# Get the residue names and residue ids of the nucleotides
residue_ids = []
residue_names = []
for residue in struc.residue_iter(nucleotides):
    mapped_nucleotide, exact_match = struc.map_nucleotide(residue)
    if mapped_nucleotide is None:
        continue
    residue_ids.append(residue[0].res_id)
    if exact_match:
        residue_names.append(mapped_nucleotide)
    else:
        residue_names.append(mapped_nucleotide.lower())

# Compute the basepairs
base_pairs = struc.base_pairs(nucleotides)
glycosidic_bonds = struc.base_pairs_glycosidic_bond(nucleotides, base_pairs)
edges = struc.base_pairs_edge(nucleotides, base_pairs)
base_pairs = struc.get_residue_positions(
    nucleotides, base_pairs.flatten()
).reshape(base_pairs.shape)

# Get the one-letter-codes of the bases
base_labels = []
for base in struc.residue_iter(nucleotides):
    base_labels.append(base.res_name[0])

# Use the base labels to indicate the Leontis-Westhof nomenclature
annotations = []
for bases, edge_types, orientation in zip(base_pairs, edges, glycosidic_bonds):
    for base, edge in zip(bases, edge_types):
        if orientation == 1:
            annotation = "c"
        else:
            annotation = "t"
        if edge == 1:
            annotation += "W"
        elif edge == 2:
            annotation += "H"
        else:
            annotation += "S"
        annotations.append(annotation)

# Print the base pairings and the Leontis-Westhof classification
for i in range(base_pairs.shape[0]):
    edge1 = annotations[i]
    edge2 = annotations[i]
    print(
        f"{base_pairs[i, 0]} {residue_names[base_pairs[i, 0]]}"
        f" {glycosidic_bonds[i]} {edge1} - "
        f"{base_pairs[i, 1]} {residue_names[base_pairs[i, 1]]}"
        f" {glycosidic_bonds[i]} {edge2}"
    )