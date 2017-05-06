#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.stats import spearmanr, pearsonr, ttest_ind
from pandas import DataFrame, read_csv, pivot_table

# -- CNV interactions
cint = read_csv('./tables/ppairs_cnv_regulation_all.csv')
cint_set = {(px, py) for px, py in cint[['px', 'py']].values}

# -- COG complex matrix
cog = read_csv('./tables/complex_ko_cog_matrix.csv', index_col=0).T.unstack().reset_index()
cog.columns = ['px', 'py', 'wb']
cog['wb'] /= 100

# -- EIF3 complex KD
eif = read_csv('./tables/complex_ko_eif3.csv', index_col=0).unstack().reset_index()
eif.columns = ['px', 'py', 'wb']

# --
complex_ko = cog.append(eif)
complex_ko['beta'] = [cint.loc[(cint['px'] == px) & (cint['py'] == py), 'beta'].values[0] if (px, py) in cint_set else np.nan for px, py in complex_ko[['px', 'py']].values]
complex_ko['fdr'] = [cint.loc[(cint['px'] == px) & (cint['py'] == py), 'fdr'].values[0] if (px, py) in cint_set else np.nan for px, py in complex_ko[['px', 'py']].values]
complex_ko = complex_ko.dropna().sort_values('fdr')
complex_ko['signif'] = ['FDR < 5%' if i < 0.05 else 'All' for i in complex_ko['fdr']]
complex_ko.to_csv('./tables/complex_ko_table.csv', index=False)


# -- Plot
plot_df = complex_ko.dropna()

# - Scatter
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'lines.linewidth': .75})
g = sns.jointplot(
    'wb', 'beta', plot_df, 'reg', color='#34495e', space=0,
    annot_kws={'template': 'Spearman: {val:.2g}, p-value: {p:.1e}'},
    joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5, 'color': '#34495e'}},
    marginal_kws={'hist': True, 'rug': False, 'kde': False, 'bins': 15},
    # xlim=[0, 150],
    # ylim=[-0.2, 0.6],
    stat_func=spearmanr
)
g.set_axis_labels('Measured KO impact', 'Estimated KO impact')
plt.gcf().set_size_inches(3, 3)
plt.savefig('./reports/wb_complex_ko_scatter.pdf', bbox_inches='tight')
plt.close('all')

# - Swarmplot
order = ['All', 'FDR < 5%']
pal = dict(zip(*(order, ['#dfdfdf', '#34495e'])))

t, pval = ttest_ind(plot_df.loc[plot_df['signif'] == 'All', 'wb'], plot_df.loc[plot_df['signif'] == 'FDR < 5%', 'wb'], equal_var=False)

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'lines.linewidth': .75})
sns.stripplot('signif', 'wb', data=plot_df, palette=pal, split=True, edgecolor='white', lw=.3, size=3, order=order, jitter=.2)
sns.boxplot('signif', 'wb', data=plot_df, palette=pal, sym='', order=order)
plt.xlabel('')
plt.ylabel('WB quantification')
plt.title('Welch\'s t-test\np-value=%.2e' % pval)
plt.gcf().set_size_inches(1, 3)
plt.savefig('./reports/wb_complex_ko_swarmplot.pdf', bbox_inches='tight')
plt.close('all')
