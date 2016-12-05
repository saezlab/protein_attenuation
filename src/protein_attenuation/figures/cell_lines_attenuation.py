#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from protein_attenuation import palette
from scipy.stats.stats import ttest_ind
from pandas import DataFrame, Series, read_csv, pivot_table


# -- Cell lines: protoemics
brca_proteomics = read_csv('./data/brca_cell_lines_proteomics_preprocessed.csv', index_col=0)
hgsc_proteomics = read_csv('./data/hgsc_cell_lines_proteomics_preprocessed.csv', index_col=0)

# -- Cell lines: transcriptomics
transcriptomics = read_csv('./data/sanger_gene_experssion_rma.tsv', sep='\t')
transcriptomics = pivot_table(transcriptomics, index='GENE_NAME', columns='SAMPLE_NAME', values='Z_SCORE', fill_value=np.nan, aggfunc=np.mean)

# -- Cell lines: Copy-number
cnv = read_csv('./data/sanger_copy_number.tsv', sep='\t')
cnv['value'] = [1 if i == 'gain' else (-1 if i == 'low' else 0) for i in cnv['MUT_TYPE']]
cnv['gene'] = [i.split('_')[0] for i in cnv['gene_name']]

cnv = cnv.groupby(['SAMPLE_NAME', 'gene'])['value'].agg(lambda x: np.nan if len(set(x)) > 1 else list(x)[0])
cnv = cnv.reset_index()
cnv = cnv.dropna()
cnv = pivot_table(cnv, index='gene', columns='SAMPLE_NAME', values='value', fill_value=0)


# -- Overlap
samples = set(brca_proteomics).union(hgsc_proteomics).intersection(transcriptomics).intersection(cnv)
print 'Overlap: ', len(samples)


# -- Attenuated proteins in tumours
t_cors = read_csv('./tables/proteins_correlations.csv', index_col=0)
p_attenuated = set(t_cors[t_cors['cluster'] == 1].index)

cors = []
for g in set(brca_proteomics.index).intersection(transcriptomics.index).intersection(cnv.index):
    df = DataFrame({'CNV': cnv.ix[g, samples], 'Transcriptomics': transcriptomics.ix[g, samples], 'Proteomics': brca_proteomics.ix[g, samples]}).dropna().corr()

    cors.append({
        'gene': g,
        'cnv_tran': df.ix['CNV', 'Transcriptomics'],
        'cnv_prot': df.ix['CNV', 'Proteomics'],
        'type': 'BRCA'
    })

for g in set(hgsc_proteomics.index).intersection(transcriptomics.index).intersection(cnv.index):
    df = DataFrame({'CNV': cnv.ix[g, samples], 'Transcriptomics': transcriptomics.ix[g, samples], 'Proteomics': hgsc_proteomics.ix[g, samples]}).dropna().corr()

    cors.append({
        'gene': g,
        'cnv_tran': df.ix['CNV', 'Transcriptomics'],
        'cnv_prot': df.ix['CNV', 'Proteomics'],
        'type': 'HGSC'
    })

cors = DataFrame(cors).dropna()
cors['diff'] = cors['cnv_tran'] - cors['cnv_prot']
cors['tumour'] = [
    'Not attenuated' if i not in p_attenuated else ('Attenuated (> %.1f)' % (np.floor((t_cors.ix[i, 'diff'] if t_cors.ix[i, 'diff'] < .5 else .5) * 10) / 10)) for i in cors['gene']
]
print cors.sort('diff')


# -- Plot
t, pval = ttest_ind(cors[cors['tumour'] != 'Not attenuated']['diff'], cors[cors['tumour'] == 'Not attenuated']['diff'])

order = ['Not attenuated', 'Attenuated (> 0.2)', 'Attenuated (> 0.3)', 'Attenuated (> 0.4)', 'Attenuated (> 0.5)']
pal = dict(zip(*(order, ['#99A3A4'] + sns.light_palette('#E74C3C', 5).as_hex()[1:])))

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(cors, size=1.25, aspect=2.5)
g = g.map_dataframe(sns.stripplot, 'diff', 'type', 'tumour', orient='h', size=2, jitter=.2, alpha=.4, linewidth=.1, edgecolor='white', palette=pal, hue_order=order, split=True)
g = g.map_dataframe(sns.boxplot, 'diff', 'type', 'tumour', orient='h', linewidth=.3, sym='', palette=pal, hue_order=order)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Attenuation score (cell lines)', '')
g.add_legend(title='Tumour attenuation')
plt.title('Protein attenuation validation\np-value: %.2e' % pval)
g.despine(trim=True)
plt.savefig('./reports/protein_attenuation_validation_boxplots.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/protein_attenuation_validation_boxplots.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
