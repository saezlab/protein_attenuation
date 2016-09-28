#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves
#
# Heatmap annotation thanks to Marco Galardini:
# https://github.com/mgalardini/python_plotting_snippets/blob/master/notebooks/clusters.ipynb

import numpy as np
import seaborn as sns
import itertools as it
import matplotlib.pyplot as plt
import fastcluster as fst
import matplotlib.patches as patches
from scipy.cluster import hierarchy
from scipy.stats.mstats import mquantiles
from cptac.utils import jaccard
from pandas import DataFrame, Series, read_csv
from seaborn.external.husl import hex_to_husl
from pymist.utils.corumdb import get_complexes_pairs, get_complexes_dict, get_complexes_name
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import
# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape


# -- CORUM
uniprot = read_uniprot_genename()

corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for px, py in corum for s, t in [(px, py), (py, px)] if s in uniprot and t in uniprot}
corum = {(px, py) for px, py in corum if px in proteomics.index and py in proteomics.index}
print 'corum', len(corum)

# -- Overlap
proteins = set(proteomics.index)
samples = set(proteomics)
print 'proteins, samples', len(proteins), len(samples)

# -- CORUM dict
corum_dict = get_complexes_dict()
corum_dict = {k: {uniprot[p][0] for p in corum_dict[k] if p in uniprot} for k in corum_dict}
corum_dict = {k: corum_dict[k].intersection(proteins) for k in corum_dict if len(corum_dict[k].intersection(proteins)) > 1}

corum_n = get_complexes_name()
print 'corum_dict', len(corum_dict)

# -- Protein-pairs correlation
p_corr = proteomics.T.corr(method='pearson')
p_corr.columns.name = 'Proteins'
p_corr.index.name = 'Proteins'
print 'p_corr', p_corr.shape


# -- Protein correlation heatmap
cmap = sns.diverging_palette(h_neg=hex_to_husl('#00B4FE')[0], h_pos=hex_to_husl('#FE4A00')[0], as_cmap=True, center='light', sep=10, l=50)

sns.set(style='white', font_scale=.75)
sns.clustermap(p_corr, figsize=(5, 5), cmap=cmap, center=0, vmax=1, vmin=-1, xticklabels=False, yticklabels=False)
plt.suptitle('Protein pairwise correlation')
plt.savefig('./reports/protein_clustering_heatmap.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'

#
sns.set(style='white', font_scale=1.)
cmap = sns.diverging_palette(h_neg=hex_to_husl('#00B4FE')[0], h_pos=hex_to_husl('#FE4A00')[0], as_cmap=True, center='light', sep=10, l=50)

for c in [387, 3544]:
    plot_df = p_corr.ix[corum_dict[c], corum_dict[c]]

    sns.clustermap(plot_df, figsize=(5, 5), cmap=cmap, center=0, vmax=1, vmin=-1, annot=True, fmt='.2f')
    plt.suptitle(corum_n[c])
    plt.savefig('./reports/protein_clustering_heatmap_%d.png' % c, bbox_inches='tight', dpi=300)
    plt.close('all')
print '[INFO] Done'


# -- Clustering analysis
lkg = fst.linkage(p_corr, method='average')

n_steps = 50
t_min, t_max = min(lkg[:, 2]), 0.7 * max(lkg[:, 2])
t_step = (t_max - t_min) / n_steps
print 'Thresholds (%d): min %.3f, max %.3f, step %.3f' % (n_steps, t_min, t_max, t_step)

thresolds = {}
for thres in np.arange(t_min, t_max, t_step):
    clusters = hierarchy.fcluster(lkg, thres, 'distance')
    clusters = {c: set(p_corr.index[clusters == c]) for c in set(clusters)}
    clusters = {k: clusters[k] for k in clusters if len(clusters[k]) > 1}

    pred = {(px, py) for k in clusters for px, py in it.permutations(clusters[k], 2)}

    # Jaccard
    jaccard_score = jaccard(corum, pred)

    # F1
    tp = len(corum.intersection(pred))
    fn = len(corum.difference(pred))
    fp = len(pred.difference(corum))

    precision = float(tp) / (tp + fp)
    recall = float(tp) / (tp + fn)

    f1 = 2 * ((precision * recall) / (precision + recall)) if precision != 0 and recall != 0 else 0

    # Store
    thresolds[thres] = {'jaccard': jaccard_score, 'f1': f1, 'precision': precision, 'recall': recall}

thresolds = DataFrame(thresolds).T
print thresolds.sort('f1', ascending=False)


# -- Clusters vs Complexes similarity
thres = thresolds.argmax()

clusters = hierarchy.fcluster(lkg, thres, 'distance')
clusters = {c: set(p_corr.index[clusters == c]) for c in set(clusters)}
clusters = {k: clusters[k] for k in clusters if len(clusters[k]) > 1}

similarity = [{'cluster': cluster, 'complex': p_complex, 'jaccard': jaccard(clusters[cluster], corum_dict[p_complex])} for cluster in clusters for p_complex in corum_dict]
similarity = DataFrame(similarity)
similarity['name'] = [corum_n[c] for c in similarity['complex']]
print similarity.sort('jaccard', ascending=False)

# _ subset
similarity_larger = similarity[[len(clusters[c1]) > len(corum_dict[c2]) for c1, c2 in similarity[['cluster', 'complex']].values]]
print similarity_larger.sort('jaccard', ascending=False)

cluster, p_complex = 1229, 162
cluster, p_complex = clusters[cluster], corum_dict[p_complex]
print 'Cluster: %s' % '; '.join(cluster)
print 'Complex: %s' % '; '.join(p_complex)
print 'Difference: %s' % '; '.join(cluster.difference(p_complex))

# Heatmap
cmap = sns.diverging_palette(h_neg=hex_to_husl('#00B4FE')[0], h_pos=hex_to_husl('#FE4A00')[0], as_cmap=True, center='dark', sep=20, l=50)

sns.set(style='white', font_scale=.75)
sns.clustermap(p_corr.ix[cluster, cluster], figsize=(5, 5), cmap=cmap, center=0, vmax=1, vmin=-1, linewidths=.1, fmt='.2f', annot=True)
plt.savefig('./reports/protein_clustering_heatmap_example.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# # -- Heatmap
# thres = thresolds.argmax()
#
# cmap = sns.diverging_palette(h_neg=hex_to_husl('#00B4FE')[0], h_pos=hex_to_husl('#FE4A00')[0], as_cmap=True, center='dark', sep=20, l=50)
#
# sns.set(style='white', font_scale=.75)
# mclust = sns.clustermap(p_corr, xticklabels=False, yticklabels=False, figsize=(10, 10), cmap=cmap, center=0, vmax=1, vmin=-1, row_linkage=lkg, col_linkage=lkg)
# mclust.ax_col_dendrogram.hlines(
#     thres,
#     0,
#     p_corr.shape[0]*10,
#     colors='g',
#     linewidths=2,
#     zorder=1
# )
#
# mclust.ax_row_dendrogram.vlines(
#     thres,
#     0,
#     p_corr.shape[0]*10,
#     colors='g',
#     linewidths=2,
#     zorder=1
# )
#
# # Extract the clusters
# clusters = hierarchy.fcluster(lkg, thres, 'distance')
# for c in set(clusters):
#     # Retrieve the position in the clustered matrix
#     index = [x for x in range(p_corr.shape[0]) if mclust.data2d.columns[x] in p_corr.index[clusters == c]]
#
#     # No singletons, please
#     if len(index) == 1:
#         continue
#
#     # Draw a rectangle around the cluster
#     mclust.ax_heatmap.add_patch(
#         patches.Rectangle(
#             (min(index),
#              p_corr.shape[0] - max(index) - 1),
#                 len(index),
#                 len(index),
#                 facecolor='none',
#                 edgecolor='g',
#                 lw=1.5)
#         )
#
# plt.suptitle('Cluster matrix')
#
# plt.savefig('./reports/protein_clustering_heatmap.png', bbox_inches='tight', dpi=150)
# plt.close('all')
# print '[INFO] Done'
