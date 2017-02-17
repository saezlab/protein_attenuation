#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from pandas import DataFrame, Series, read_csv


# -- Import
phoshoproteomics = read_csv('./data/cptac_phosphoproteomics_corrected_normalised.csv', index_col=0)
proteomics = read_csv('./data/cptac_proteomics_tmt_corrected_normalised.csv', index_col=0)


# -- Overlap
samples = set(proteomics).intersection(phoshoproteomics)
proteins = set(proteomics.index)
psites = {p for p in phoshoproteomics.index if p.split('_')[0] in proteins}


# -- Correlations
# p = 'EPS8L1_T187'
def phospho_correlation(p):
    print p
    df = DataFrame({'phospho': phoshoproteomics.ix[p], 'prot': proteomics.ix[p.split('_')[0]]}).dropna()
    return pearsonr(df['prot'], df['phospho'])[0] if df.shape[0] > (len(samples) * .5) else np.nan

correlation = Series({p: phospho_correlation(p) for p in psites}).dropna()
print correlation.sort_values()


# -- Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.distplot(correlation, bins=30, color='#bfbfbf', kde_kws={'lw': .3}, hist_kws={'lw': .3})
plt.axvline(0, ls='--', lw=0.1, c='black', alpha=.3)
sns.despine(trim=True)
plt.ylabel('Density')
plt.xlabel('Pearson')
plt.title('Phospho-site ~ protein correlation\n(mean: %.2f)' % correlation.mean())
plt.gcf().set_size_inches(3, 2)
plt.savefig('./reports/psites_protein_correlation_histogram.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/psites_protein_correlation_histogram.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] PCA after normalisation: ', './reports/psites_protein_correlation_histogram.pdf'
