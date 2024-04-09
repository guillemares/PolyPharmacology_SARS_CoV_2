#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import linkage, cophenet, fcluster, dendrogram
from superimposition import superimposition

def heatmap(data_rmsd, data_tm):
    """
    Given the RMSD and TM scores dataframes
    from superimposition.py, create a heatmap.
    """
    sns.set_theme()
    plt.figure(figsize=(18, 18))
    sns.heatmap(data_rmsd, annot=True, fmt=".2f", annot=True, cmap="YlGnBu", fmt=".2f", linewidths=0.5, linecolor="white")
    plt.title("Superimposition Heatmap", fontsize=22)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()
    # plt.savefig("heatmap.pdf")

def clustermap(data_rmsd, data_tm):
    """
    Given the RMSD and TM scores dataframes
    from superimposition.py, create a clustermap.
    """
    rmsd_cluster = sns.clustermap(data_rmsd, cmap="YlGnBu", linewidths=0.5, linecolor="white", figsize=(12, 12), cbar_pos=(0.075, 0.24, 0.025, 0.6), row_cluster=True, col_cluster=True)
    plt.title("Superimposition Clustermap", fontsize=22)
    for dendrogram in rmsd_cluster.ax_row_dendrogram.collections:
        dendrogram.set_visible(False)
    for dendrogram in rmsd_cluster.ax_col_dendrogram.collections:
        dendrogram.set_linewidth(1.2)
    plt.show()
    # plt.savefig("clustermap.pdf")

def plot_dendrograms(self, data_rmsd, data_tm):
    """
    Given the RMSD and TM scores dataframes
    from superimposition.py, create dendrograms.
    """
    plt.figure(figsize=(25, 15))
    plt.title('Hierarchical Clustering Dendrogram', fontsize=20)
    plt.xlabel('sample index', fontsize=18)
    plt.ylabel('distance', fontsize=18)
    dendrogram = sch.dendrogram(sch.linkage(data_rmsd.values, method='average', metric='euclidean'), labels=data_rmsd.index, leaf_rotation=50., leaf_font_size=8.)
    plt.yticks(fontsize=10)
    plt.savefig("dendrogram.png", dpi=300)

    return 0

def clustering(self, data_rmsd, data_tm):
    """
    Given the RMSD and TM scores dataframes
    from superimposition.py, perform clustering.
    """
    clusters = fcluster(sch.linkage(data_rmsd.values, method='average', metric='euclidean'), 0.8, criterion='distance')
    clusters_dic = {}
    for i, cluster in enumerate(clusters):
        if cluster not in clusters_dic:
            clusters_dic[cluster] = []
        clusters_dic[cluster].append(data_rmsd.index[i])

    for key in clusters_dic.keys():
        print('cluster%d: %s' % (key, clusters_dic[key]))

    cluster_centers = {}
    for cluster in clusters_dic.keys():
        cluster_centers[cluster] = np.mean(data_rmsd.loc[clusters_dic[cluster], clusters_dic[cluster]].values)
        print('center cluster%d: %.2f' % (cluster, cluster_centers[cluster]))

    self.clusters_dic = clusters_dic
    self.clusters = clusters
    self.cluster_centers = cluster_centers

    return 0

def closest_to_center(self, data_rmsd, data_tm):
    """
    Given the RMSD and TM scores dataframes
    from superimposition.py, find the closest
    structure to the center of the cluster.
    """

    for cluster in self.clusters_dic.keys():
        print('cluster%d: %s' % (cluster, self.clusters_dic[cluster]))
        center=self.cluster_centers[cluster]
        elements=self.clusters_dic[cluster]
        submatrix=self.data_tm_score.loc[elements, elements]
        distances=np.sum(submatrix, axis=0)
        index=np.argmin(distances)
        print('representative element: %s' % elements[index])

    return 0

def main():
    superimpobj = superimposition(method='1vs1', dir='../test/superimposition/', test=True)
    superimpobj._1_vs_1()
    heatmap(superimpobj.data_rmsd, superimpobj.data_tm_score)
    clustermap(superimpobj.data_rmsd, superimpobj.data_tm_score)
    plot_dendrograms(superimpobj.data_rmsd, superimpobj.data_tm_score)
    clustering(superimpobj.data_rmsd, superimpobj.data_tm_score)
    closest_to_center(superimpobj.data_rmsd, superimpobj.data_tm_score)

main()
