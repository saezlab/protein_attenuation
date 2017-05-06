#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pandas import Series, DataFrame, read_csv, pivot_table

# -- CNV
cnv = read_csv('./data/ccle_cna_mat.txt', sep='\t', index_col=0)

# -- RPPA
cmap = read_csv('./data/ccle_shared_cna_samples.txt', sep='\t', index_col=0)['ccle']

rppa = read_csv('./data/MCLP-v1.1-Level4.tsv', sep='\t', index_col=0).T
rppa = rppa.loc[:, [c in cmap.index for c in rppa]]
rppa.columns = cmap.ix[rppa.columns].values

# -- Transcriptomics
cmap = read_csv('./data/CCLE_GDSC_cellLineMappoing.csv', index_col=1)['CCLE name']

transcriptomics = read_csv('./data/sanger_gene_experssion_rma_preprocessed.csv', index_col=0)
transcriptomics = transcriptomics.loc[:, [c in cmap.index for c in transcriptomics]]
transcriptomics.columns = cmap.ix[transcriptomics.columns].values

# -- Overlap
samples, proteins = set(rppa).intersection(cnv).intersection(transcriptomics), set(rppa.index).intersection(cnv.index).intersection(transcriptomics.index)
print(len(samples), len(proteins))

# -- Attenuated proteins in tumours
t_cors = read_csv('./tables/protein_attenuation_table.csv', index_col=0)
p_attenuated = set(t_cors[t_cors['attenuation_potential'] == 'High'].index)

# -- Correlations
cors = []
for g in proteins:
    df = DataFrame({'CNV': cnv.ix[g, samples], 'Transcriptomics': transcriptomics.ix[g, samples], 'RPPA': rppa.ix[g, samples]}).dropna().corr()

    cors.append({
        'gene': g,
        'cnv_tran': df.ix['CNV', 'Transcriptomics'],
        'cnv_rppa': df.ix['CNV', 'RPPA'],
        'type': 'RPPA'
    })

cors = DataFrame(cors).dropna()
cors['attenuation'] = cors['cnv_tran'] - cors['cnv_rppa']
cors['attenuation_potential'] = [
    'Low' if i not in p_attenuated else ('High (> %.1f)' % (np.floor((t_cors.ix[i, 'attenuation'] if t_cors.ix[i, 'attenuation'] < .5 else .5) * 10) / 10)) for i in cors['gene']
]

# -- Plot
order = ['Low', 'High (> 0.2)', 'High (> 0.3)', 'High (> 0.4)', 'High (> 0.5)']
pal = dict(zip(*(order, ['#99A3A4'] + sns.light_palette('#E74C3C', 5).as_hex()[1:])))

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(cors, size=1.25, aspect=2.5, legend_out=True)
g = g.map_dataframe(sns.stripplot, 'attenuation', 'type', 'attenuation_potential', orient='h', size=2, jitter=.2, alpha=.4, linewidth=.1, edgecolor='white', palette=pal, hue_order=order, split=True)
g = g.map_dataframe(sns.boxplot, 'attenuation', 'type', 'attenuation_potential', orient='h', linewidth=.3, sym='', palette=pal, hue_order=order)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.despine(trim=True)
g.set_axis_labels('Protein attenuation (cell lines)', '')
g.add_legend(title='Protein attenuation\n(tumours)')
plt.setp(g._legend.get_title(), multialignment='center')
plt.title('Protein attenuation transfer from tumours to cell lines')
plt.savefig('./reports/rppa_attenuation_validation_boxplots.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/rppa_attenuation_validation_boxplots.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Protein attenuation potential cell lines: ', './reports/rppa_attenuation_validation_boxplots.pdf'
