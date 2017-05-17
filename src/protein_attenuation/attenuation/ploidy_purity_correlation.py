#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from pandas import DataFrame, Series, read_csv

# -- Imports
# - Ploidy and purity
pp = read_csv('./tables/ascat_acf_ploidy.tsv', sep='\t')

# Parse ids
pp['Sample'] = [i[:12].replace('.', '-') for i in pp['Sample']]

# Discard duplicated samples
d_samples = Series(dict(zip(*(np.unique(pp['Sample'], return_counts=True))))).sort_values()
d_samples = list(d_samples[d_samples != 1].index)

pp = pp[[i not in d_samples for i in pp['Sample']]].set_index('Sample')

# - Import attenuation
attenuation = read_csv('./tables/sample_attenuation_table.csv', index_col=0)

# -- Overlap
samples = set(pp.index).intersection(attenuation.index)

# -- Plot: ploidy
plot_df = DataFrame({
    'ploidy': pp.ix[samples, 'Ploidy'], 'attenuation': attenuation.ix[samples, 'attenuation']
})

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(
    'ploidy', 'attenuation', plot_df, 'reg', color='#34495e', space=0,
    marginal_kws={'hist': False, 'rug': False},
    annot_kws={'template': 'Pearson: {val:.2g}, p-value: {p:.1e}', 'loc': 4},
    joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5, 'alpha': .5}},
    stat_func=pearsonr
)
plt.axhline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Ploidy estimate', 'Attenuation potential')
plt.gcf().set_size_inches(3, 3)
plt.savefig('./reports/ploidy_attenuation_corrplot.pdf', bbox_inches='tight')
plt.close('all')
print('[INFO] Corr plotted!')

# -- Plot: purity
plot_df = DataFrame({
    'purity': pp.ix[samples, 'Aberrant_Cell_Fraction(Purity)'], 'attenuation': attenuation.ix[samples, 'attenuation']
})

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(
    'purity', 'attenuation', plot_df, 'reg', color='#34495e', space=0,
    marginal_kws={'hist': False, 'rug': False},
    annot_kws={'template': 'Pearson: {val:.2g}, p-value: {p:.1e}', 'loc': 4},
    joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5, 'alpha': .5}},
    stat_func=pearsonr
)
plt.axhline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Purity estimate', 'Attenuation potential')
plt.gcf().set_size_inches(3, 3)
plt.savefig('./reports/purity_attenuation_corrplot.pdf', bbox_inches='tight')
plt.close('all')
print('[INFO] Corr plotted!')
