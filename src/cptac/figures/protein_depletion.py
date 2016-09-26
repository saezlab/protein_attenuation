#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.stats import ttest_ind
from pandas import read_csv, DataFrame
from pymist.enrichment.gsea import gsea
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)
print 'cnv', cnv.shape

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape


# -- Overlap
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
genes = set(cnv.index).intersection(transcriptomics.index).intersection(proteomics.index)
print 'samples', 'genes', len(samples), len(genes)


# -- Correlations
res = {}
for g in genes:
    df = DataFrame({'CNV': cnv.ix[g, samples], 'Transcriptomics': transcriptomics.ix[g, samples], 'Proteomics': proteomics.ix[g, samples]}).dropna().corr()

    res[g] = {
        'CNV_Transcriptomics': df.ix['CNV', 'Transcriptomics'],
        'CNV_Proteomics': df.ix['CNV', 'Proteomics'],
        'Transcriptomics_Proteomics': df.ix['Transcriptomics', 'Proteomics']
    }

    print g

res = DataFrame(res).T
print res


# -- Enrichment
uniprot = read_uniprot_genename()

corum = get_complexes_pairs()
corum = {uniprot[g][0] for p in corum for g in p if g in uniprot}.intersection(genes)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = set(ppairs_cnv[ppairs_cnv['fdr'] < .05]['px'])

res['diff'] = res['CNV_Transcriptomics'] - res['CNV_Proteomics']
res['Interaction'] = ['Complex' if i in corum else 'All' for i in res.index]

# ttest_ind(res.loc[res['Interaction'] == 'Complex', 'diff'], res.loc[res['Interaction'] == 'All', 'diff'], equal_var=False)

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
sns.boxplot(x='Interaction', y='diff', data=res, linewidth=1., notch=True, sym='')
sns.stripplot(x='Interaction', y='diff', data=res, linewidth=1., alpha=.5, jitter=True, edgecolor='white')
sns.despine(trim=True)
plt.ylabel('Pearson(CNV ~ Transcriptomics) - Pearson(CNV ~ Proteomics)')
plt.gcf().set_size_inches(3, 5)
# plt.savefig('./reports/correlation_difference_violinplots.pdf', bbox_inches='tight')
plt.savefig('./reports/correlation_difference_violinplots.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Plot done'
