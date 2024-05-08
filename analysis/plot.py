import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import glob
# sys.path.insert(0,'../')
from rna import RNA

"""
Given an ensemble of RNA structures with the same
sequence, plot the interactions between the residues
previously defined in the RNA object.
"""

def plot_ensemble_interactions(dataframe, rnas):
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
    # Change the residue number to another index in the column
    # BaseId1 and BaseId2 of the dataframe
    # merged_df = init_index(merged_df)
    merged_df = dataframe
    merged_df.to_csv('merged_df.txt', index=True, sep='\t')
    for index, row in merged_df.iterrows():
        bp_label = f"{row['BaseName1']}{row['BaseId1']}-{row['BaseName2']}{row['BaseId2']}"
        if bp_label not in pair_df_index:
            pair_df_index[bp_label] = len(pair_df)
            pair_df.append(bp_label)
    def custom_sort(x):
        return int(x.split('-')[0][1:])
    pair_df = sorted(pair_df, key=custom_sort)
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(38, 25))
    x = np.arange(len(pair_df))
    width = 0.5
    unique_interactions = merged_df['Interaction Type'].unique()
    total_interactions_per_pair = np.zeros((len(pair_df), len(unique_interactions)))
    for i, pair in enumerate(pair_df):
        total_interactions = []
        for j, interaction in enumerate(unique_interactions):
            count = merged_df[(merged_df['BaseName1'] + merged_df['BaseId1'].astype(str) + "-" + merged_df['BaseName2'] + merged_df['BaseId2'].astype(str) == pair) &
                               (merged_df['Interaction Type'] == interaction)].shape[0]
            total_interactions_per_pair[i, j] = count
    bottom = None
    for i, interaction in enumerate(unique_interactions):
        color = color_palette.get(interaction, 'black')
        plt.bar(x, total_interactions_per_pair[:, i], width, label=interaction, bottom=bottom, color=color)
        if bottom is None:
            bottom = total_interactions_per_pair[:, i]
        else:
            bottom += total_interactions_per_pair[:, i]
    plot_name = f'interaction_barplot.png'
    plt.xlabel('Base Pair', fontsize=26)
    plt.ylabel('Count', fontsize=26)
    plt.title('Number of interactions per base pair', fontsize=28)
    plt.xticks(x, pair_df, fontsize=18, rotation=45)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.savefig(plot_name)
    plt.show()

    return 0

def color_palette():
    """
    Create a dictionary with the interaction types and the color
    assigned to each interaction type.
    The colors are chosen manually to be visually distinguishable
    and easy to change when needed.
    """
    possible_interactions = {

        'cWcW': 'cornflowerblue', 'cWcH': 'royalblue', 'cWcS': 'blue', 'cWcX': 'darkblue',
        'cHcW': 'green', 'cHcH': 'limegreen', 'cHcS': 'springgreen', 'cHcX': 'darkgreen',
        'cScW': 'slateblue', 'cScH': 'rebeccapurple', 'cScS': 'mediumorchid', 'cScX': 'violet',
        'cXcW': 'deeppink', 'cXcH': 'pink', 'cXcS': 'hotpink', 'cXcX': 'salmon',

        'tWtW': 'aquamarine', 'tWtH': 'mediumturquoise', 'tWtS': 'paleturquoise', 'tWtX': 'cyan',
        'tHtW': 'steelblue', 'tHtH': 'lightslategray', 'tHtS': 'grey', 'tHtX': 'silver',
        'tStW': 'bisque', 'tStH': 'burlywood', 'tStS': 'antiquewhite', 'tStX': 'blanchedalmond',
        'tXtW': 'goldenrod', 'tXtH': 'gold', 'tXtS': 'khaki', 'tXtX': 'darkkhaki',

        'xWxW': 'olive', 'xWxH': 'yellow', 'xWxS': 'yellowgreen', 'xWxX': 'lightyellow',
        'xHxW': 'palegreen', 'xHxH': 'darkseagreen', 'xHxS': 'chartreuse', 'xHxX': 'mediumspringgreen',
        'xSxW': 'peru', 'xSxH': 'peachpuff', 'xSxS': 'sienna', 'xSxX': 'saddlebrown',
        'xXxW': 'orangered', 'xXxH': 'red', 'xXxS': 'tomato', 'xXxX': 'crimson',

        'WobbleGU': 'firebrick'

    }

    return possible_interactions


if __name__ == "__main__":
    rnas_files = glob.glob('*.pdb')
    rnas = []
    merged_df = pd.DataFrame()
    for rna_file in rnas_files:
        rnaobj = RNA(pdbfile=rna_file)
        rnaobj.get_full_df()
        merged_df = rnaobj.merge_df(merged_df)
        print(merged_df)

    color_palette = color_palette()
    plot_ensemble_interactions(dataframe=merged_df, rnas=rnas)