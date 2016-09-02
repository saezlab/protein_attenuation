#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import palette, default_color
from matplotlib.gridspec import GridSpec
from sklearn.metrics.ranking import roc_curve, auc
from pandas import DataFrame, Series, read_csv, concat


# -- Import
pslist = read_csv('./tables/PN-ontology_MB_082016.csv').dropna(how='all', axis=1)
pslist_dict = pslist.groupby('Concatenation')['Gene'].agg(lambda x: set(x)).to_dict()
print 'pslist', pslist.shape

proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape

transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape

# -- Overlap
proteins, samples = set(pslist['Gene']).intersection(proteomics.index).intersection(transcriptomics.index), set(proteomics).intersection(transcriptomics)
print 'proteins', 'samples', len(proteins), len(samples)


# -- Protein-protein correlation
cor_df = {}
# name, d_df = 'BRCA', brca
for name, d_df in [('Proteomics', proteomics), ('Transcriptomics', transcriptomics)]:
    print name

    # -- Generate p-p correlations
    df = d_df.ix[proteins, samples].T.corr(method='pearson')
    df.values[np.tril_indices(df.shape[0], 0)] = np.nan
    df.index.name = None
    df = df.unstack().reset_index().dropna()
    df.columns = ['p1', 'p2', 'cor']

    # -- Annotate p-p interactions
    for k in pslist_dict:
        df[k] = [1 if p1 in pslist_dict[k] and p2 in pslist_dict[k] else 0 for p1, p2 in df[['p1', 'p2']].values]

    # -- Store cor
    cor_df[name] = df.copy()
print 'cor_df', len(cor_df)


# --
roc_df = []
for k in pslist_dict:
    for n in cor_df:
        curve_fpr, curve_tpr, _ = roc_curve(cor_df[n][k], cor_df[n]['cor'])
        curve_auc = auc(curve_fpr, curve_tpr)

        res = {'group': k, 'omics': n, 'fpr': curve_fpr, 'tpr': curve_tpr, 'auc': curve_auc, 'len': cor_df[n][k].sum()}

        roc_df.append(res)

roc_df = DataFrame(roc_df).dropna(subset=['auc'])
roc_df[roc_df['len'] > 1].sort('auc', ascending=False)[['group', 'omics', 'auc', 'len']].to_csv('./tables/psn_aucs.csv', index=False)
# roc_df = read_csv('./tables/psn_aucs.csv')
print roc_df.sort('auc')[['auc', 'omics', 'len', 'group']]

# Scatter
plot_df = roc_df[roc_df['len'] > 1]
plot_df = DataFrame({
    'Transcriptomics': plot_df[plot_df['omics'] == 'Transcriptomics'].set_index('group')['auc'],
    'Proteomics': plot_df[plot_df['omics'] == 'Proteomics'].set_index('group')['auc']}
).reset_index()

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})

plt.scatter(plot_df['Transcriptomics'], plot_df['Proteomics'], c=default_color)
plt.plot([0, 1], [0, 1], 'k--', lw=.3)
# plt.legend(loc='lower right')
plt.xlabel('Transcriptomics (AUC)')
plt.ylabel('Proteomics (AUC)')
plt.xlim(0, 1)
plt.ylim(0, 1)
sns.despine(trim=True, top=False, right=False)
plt.title('Proteostasis network')

plt.gcf().set_size_inches(4, 4)
plt.savefig('./reports/psn_proteomics_transcriptomics_auc_scatter.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'

# Barplot
plot_df = roc_df[roc_df['len'] > 1].sort_values('auc', ascending=False)[['group', 'auc', 'omics', 'len']].copy().reset_index()
plot_df['len'] = np.log10(plot_df['len'])

sns.set(style='white', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .0, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'}, font_scale=0.75)
fig, axes = plt.subplots(ncols=2, sharey=True)
pos = 0

for f in ['auc', 'len']:
    sns.stripplot(x=f, y='group', hue='omics', data=plot_df, palette=palette, linewidth=0, size=4, orient='h', split=False, ax=axes[pos], jitter=True)

    axes[pos].xaxis.grid(False)
    axes[pos].yaxis.grid(True)

    sns.despine(bottom=True, left=True, ax=axes[pos])

    axes[pos].legend(loc=4)

    axes[pos].set_xlabel('AUC' if f == 'auc' else 'Number of pairs associated (log10)')

    pos += 1


plt.gcf().set_size_inches(6, 8)
plt.savefig('./reports/psn_proteomics_transcriptomics_auc_barplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
