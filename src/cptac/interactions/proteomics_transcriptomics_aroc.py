#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import igraph
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import matplotlib.pyplot as plt
from cptac import wd, palette, default_color
from matplotlib.gridspec import GridSpec
from pymist.utils.stringdb import get_stringdb
from pymist.utils.biogriddb import get_biogriddb
from sklearn.metrics.ranking import roc_curve, auc
from pymist.utils.corumdb import get_complexes_pairs
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import imputed proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print 'proteomics', proteomics.shape

transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print 'transcriptomics', transcriptomics.shape


# -- Import protein-protein interactions
uniprot = read_uniprot_genename()

# CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print 'corum', len(corum)

# String
string = get_stringdb(900)
string = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in string for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print 'string', len(string)

# Biogrid
biogrid = get_biogriddb()
biogrid = {(s, t) for p1, p2 in biogrid for s, t in [(p1, p2), (p2, p1)]}
print 'biogrid', len(biogrid)

# Omnipath
omnipath = read_csv('%s/files/omnipath_interactions.txt' % wd, sep='\t')
omnipath = {(uniprot[s][0], uniprot[t][0]) for s, t in omnipath[['source', 'target']].values if s in uniprot and t in uniprot}
print 'omnipath', len(omnipath)

omnipath_action = read_csv('%s/files/omnipath_interactions.txt' % wd, sep='\t')
omnipath_sources = {s for i in omnipath_action['sources'] for s in i.split(';')}
omnipath_action = {s: {(uniprot[s][0], uniprot[t][0]) for s, t in omnipath_action[[s in i for i in omnipath_action['sources']]][['source', 'target']].values if s in uniprot and t in uniprot} for s in omnipath_sources}
print 'omnipath_action', len(omnipath_action)


# -- Overlap
proteins, samples = set(proteomics.index).intersection(transcriptomics.index), set(proteomics).intersection(transcriptomics)
print 'proteins', 'samples', len(proteins), len(samples)


# -- Protein-protein correlation
dfs = []
for name, d_df in [('Proteomics', proteomics), ('Transcriptomics', transcriptomics)]:
    df = d_df.ix[proteins, samples].T.corr(method='pearson')
    df.values[np.tril_indices(df.shape[0], 0)] = np.nan
    df.index.name = None
    df = df.unstack().reset_index().dropna()
    df.columns = ['p1', 'p2', name]

    df = df.set_index(['p1', 'p2'])

    dfs.append(df)

dfs = concat(dfs, axis=1).reset_index()
dfs['Sum'] = dfs['Proteomics'] + dfs['Transcriptomics']
print dfs


# - Annotate protein-pairs
dfs['CORUM'] = [1 if ((p1, p2) in corum) else 0 for p1, p2 in dfs[['p1', 'p2']].values]
dfs['STRING'] = [1 if ((p1, p2) in string) else 0 for p1, p2 in dfs[['p1', 'p2']].values]
dfs['BioGRID'] = [1 if ((p1, p2) in biogrid) else 0 for p1, p2 in dfs[['p1', 'p2']].values]
dfs['OmniPath'] = [1 if ((p1, p2) in omnipath) else 0 for p1, p2 in dfs[['p1', 'p2']].values]

# for action in omnipath_action:
#     df['OmniPath_%s' % action] = [1 if ((p1, p2) in omnipath_action[action]) else 0 for p1, p2 in df[['p1', 'p2']].values]
#     print 'OmniPath_%s' % action
print dfs


# -- Calculate ROC curves
cor_res = {}
for i in ['Proteomics', 'Transcriptomics', 'Sum']:
    cor_res[i] = {}
    for db in ['CORUM', 'STRING', 'BioGRID', 'OmniPath']:
        curve_fpr, curve_tpr, _ = roc_curve(dfs[db], dfs[i])
        cor_res[i][db] = (curve_fpr, curve_tpr)

# for action in omnipath_action:
#     curve_fpr, curve_tpr, _ = roc_curve(df['OmniPath_%s' % action], df['cor'])
#     cor_res['OmniPath_%s' % action] = (curve_fpr, curve_tpr)
#
print 'cor_res', len(cor_res)


# -- Plot
palette['Sum'] = default_color

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
fig, gs, pos = plt.figure(figsize=(13, 3)), GridSpec(1, 4, hspace=.3, wspace=.3), 0

for db in ['CORUM', 'STRING', 'BioGRID', 'OmniPath']:
    ax = plt.subplot(gs[pos])

    for i in ['Proteomics', 'Transcriptomics', 'Sum']:
        curve_fpr, curve_tpr = cor_res[i][db]
        plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (i, auc(curve_fpr, curve_tpr)), c=palette[i])

    ax.plot([0, 1], [0, 1], 'k--', lw=.3)
    sns.despine(trim=True)
    ax.legend(loc='lower right')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title(db)

    pos += 1

plt.savefig('%s/reports/proteomics_ppi_aroc_sum.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Plot OmniPtah resources ROCs
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs, pos = GridSpec(8, 5, hspace=.45, wspace=.45), 0

for action in omnipath_action:
    ax = plt.subplot(gs[pos])

    for name in ['BRCA', 'COREAD', 'HGSC', 'Proteomics', 'Transcriptomics']:
        curve_fpr, curve_tpr = cor_res[name]['OmniPath_%s' % action]
        plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (name, auc(curve_fpr, curve_tpr)), c=palette[name])

    ax.plot([0, 1], [0, 1], 'k--', lw=.3)
    sns.despine(trim=True)
    ax.legend(loc='lower right')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title('OmniPath \nresource: %s' % action)

    pos += 1

plt.gcf().set_size_inches(18, 29)
plt.savefig('%s/reports/proteomics_ppi_aroc_omnipath_type.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
