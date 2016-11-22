#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.stats import ttest_ind
from pandas import read_csv, DataFrame
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_dict


# -- Human mouse orthologs
orthologs = read_csv('./tables/human_mouse_orthologs.txt', sep='\t')
orthologs = orthologs.groupby('MGI symbol').agg(lambda x: set(x)).to_dict()['HGNC symbol']

# -- CORUM
uniprot = read_uniprot_genename()
corum_dict = get_complexes_dict()
corum_dict = {k: {uniprot[g][0] for g in corum_dict[k] if g in uniprot} for k in corum_dict}

corum_proteins = {p for k in corum_dict for p in corum_dict[k]}
print 'corum', len(corum_proteins)

# -- Improt regression results
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans_interactions = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans_interactions for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')

px, py = set(ppairs_cnv['px']), set(ppairs_cnv['py'])


# -- Import protein degradation profiles
ned = read_csv('./tables/ned_proteins.csv')[['Gene names', 'score', 'Degradation profile', 'Validation score']]
ned = ned[[i in orthologs for i in ned['Gene names']]]
ned['human'] = [';'.join([x for x in orthologs[i] if str(x) != 'nan']) for i in ned['Gene names']]
ned['type'] = ['Px' if i in px else ('Py' if i in py else ('Complex' if i in corum_proteins else 'Other')) for i in ned['human']]


# --
plot_df = ned[['human', 'type', 'Degradation profile', 'score']].copy()
plot_df = DataFrame([
    {'protein': plot_df.ix[i, 'human'], 'type': plot_df.ix[i, c], 'group': 'Degradation' if c == 'Degradation profile' else 'Complex', 'score': plot_df.ix[i, 'score']}
    for i in plot_df.index for c in ['type', 'Degradation profile']
])

t, pval = ttest_ind(ned.ix[ned['type'] == 'Px', 'score'], ned.ix[ned['type'] == 'Py', 'score'], equal_var=False)

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(plot_df, size=1, aspect=2, row='group', sharey=False)
g = g.map_dataframe(sns.stripplot, x='score', y='type', orient='h', size=2, jitter=.2, alpha=.35, linewidth=.1, edgecolor='white', color='#99A3A4')
g = g.map_dataframe(sns.boxplot, x='score', y='type', orient='h', linewidth=.3, sym='', color='#99A3A4')
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Degradation score')
g.set_titles('{row_name}')
g.despine(trim=True)
g.set_ylabels('')
plt.savefig('./reports/ned_proteins_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/ned_proteins_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'
