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

# -- Turnover rates
turnover = read_csv('./files/proteins_turnovers_preprocessed.csv').dropna(subset=['Uniprot IDs human'])
turnover = DataFrame([{'protein': i, 'p_halflife': r, 't_halflife': t} for p, r, t in turnover[['Uniprot IDs human', 'Protein half-life average [h]', 'mRNA half-life average [h]']].values for i in p.split(';')])
turnover['protein'] = [uniprot[i][0] for i in turnover['protein']]
turnover = turnover.groupby('protein').max().dropna()

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
order = ['Px', 'Py', 'Complex', 'Other']

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(ned[['type', 'score']], size=1, aspect=2)
g = g.map_dataframe(sns.stripplot, x='score', y='type', orient='h', size=2, jitter=.2, alpha=.35, linewidth=.1, edgecolor='white', color='#99A3A4', order=order)
g = g.map_dataframe(sns.boxplot, x='score', y='type', orient='h', linewidth=.3, sym='', color='#99A3A4', order=order)
g = g.map(plt.axvline, x=0, ls='--', lw=0.1, c='black', alpha=.5)
g = g.map(plt.axvline, x=ned[ned['Degradation profile'] == 'ED']['score'].max(), ls='-', lw=0.3, c='black', alpha=.7)
g.set_axis_labels('Degradation score')
g.set_titles('{row_name}')
g.despine(trim=True)
g.set_ylabels('')
plt.savefig('./reports/ned_proteins_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/ned_proteins_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'
