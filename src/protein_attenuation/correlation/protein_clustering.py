#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves
#
# Heatmap annotation thanks to Marco Galardini:
# https://github.com/mgalardini/python_plotting_snippets/blob/master/notebooks/clusters.ipynb

import numpy as np
import seaborn as sns
import fastcluster as fst
import matplotlib.pyplot as plt
from pandas import read_csv
from scipy.cluster.hierarchy import leaves_list
from protein_attenuation.utils import get_complexes_pairs, get_complexes_dict, get_complexes_name, read_uniprot_genename


# -- Import
# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)

# -- CORUM
uniprot = read_uniprot_genename()

corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for px, py in corum for s, t in [(px, py), (py, px)] if s in uniprot and t in uniprot}
corum = {(px, py) for px, py in corum if px in proteomics.index and py in proteomics.index}

# -- Overlap
proteins = set(proteomics.index)
samples = set(proteomics)

# -- CORUM dict
corum_dict = get_complexes_dict()
corum_dict = {k: {uniprot[p][0] for p in corum_dict[k] if p in uniprot} for k in corum_dict}
corum_dict = {k: corum_dict[k].intersection(proteins) for k in corum_dict if len(corum_dict[k].intersection(proteins)) > 1}

corum_n = get_complexes_name()

# -- Protein-pairs correlation
p_corr = proteomics.T.corr(method='pearson')
p_corr.columns.name = 'Proteins'
p_corr.index.name = 'Proteins'


# -- Protein correlation heatmap
plot_df = p_corr.copy()

lkg = fst.linkage_vector(plot_df.values, method='median', metric='euclidean')
order = list(plot_df.index[leaves_list(lkg)])

plot_df = plot_df.ix[order, order]

mask = np.zeros_like(plot_df)
mask[np.triu_indices_from(mask, k=0)] = True

sns.set(style='white', font_scale=1.)
g = sns.heatmap(plot_df, cmap='YlGnBu', center=0, vmax=1, vmin=-1, square=True, mask=mask, linewidths=.0, cbar=False, xticklabels=False, yticklabels=False)

plt.gcf().set_size_inches(3, 3)

plt.suptitle('Protein pairwise correlation')
plt.savefig('./reports/protein_clustering_heatmap.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Pairwise correlation heatmap: ', './reports/protein_clustering_heatmap.png'


# -- Representative complexes
sns.set(style='white', font_scale=.75)
for c, s, t in [(387, 6, 'MCM complex'), (162, 5, 'COG complex')]:
    plot_df = p_corr.ix[corum_dict[c], corum_dict[c]]

    lkg = fst.linkage_vector(plot_df.values, method='median', metric='euclidean')
    order = list(plot_df.index[leaves_list(lkg)])

    plot_df = plot_df.ix[order, order]

    mask = np.zeros_like(plot_df)
    mask[np.triu_indices_from(mask, k=1)] = True

    plt.gcf().set_size_inches(2, 2)

    # g = sns.clustermap(plot_df, figsize=(6, 6), cmap=cmap, center=0, vmax=1, vmin=-1, annot=True, fmt='.2f')
    g = sns.heatmap(plot_df, cmap='YlGnBu', center=0, vmax=1, vmin=-1, annot=True, fmt='.2f', square=True, mask=mask, linewidths=.3, annot_kws={'size': s}, cbar=False)

    plt.suptitle(t)
    plt.savefig('./reports/protein_clustering_heatmap_%d.png' % c, bbox_inches='tight', dpi=300)
    plt.savefig('./reports/protein_clustering_heatmap_%d.pdf' % c, bbox_inches='tight')
    plt.close('all')
    print '[INFO] Representative complexes: ', './reports/protein_clustering_heatmap_%d.pdf' % c


# -- Correlation heatmap color bar
sns.set(style='ticks', font_scale=1., rc={'xtick.direction': 'in', 'ytick.direction': 'in', 'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.palplot(sns.color_palette('YlGnBu', n_colors=7))
plt.title('Pearson correlation')
plt.gcf().set_size_inches(3, .5)
plt.savefig('./reports/protein_clustering_heatmap_cbar.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap color bar: ', './reports/protein_clustering_heatmap_cbar.pdf'