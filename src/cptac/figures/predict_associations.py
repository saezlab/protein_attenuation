#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, Series, read_csv
from seaborn.external.husl import hex_to_husl
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import
# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape

# CORUM
uniprot = read_uniprot_genename()

corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print 'corum', len(corum)


# -- Protein-pairs correlation
p_corr = proteomics.T.corr(method='pearson')
print p_corr


# -- Heatmap
mask = np.zeros_like(p_corr)
mask[np.triu_indices_from(mask)] = True

cmap = sns.diverging_palette(h_neg=hex_to_husl('#00B4FE')[0], h_pos=hex_to_husl('#FE4A00')[0], as_cmap=True, center='dark', sep=20, l=50)

sns.set(style='white', font_scale=.75)
sns.clustermap(p_corr.ix[:1000, :1000], col_cluster=True, row_cluster=True, xticklabels=False, yticklabels=False, figsize=(5, 5), cmap=cmap, center=0)
plt.savefig('./reports/protein_clustering_heatmap.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'
