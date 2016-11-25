#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame
from scipy.stats.stats import ttest_ind


# -- Improt regression results
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans_interactions = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans_interactions for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')

px, py = set(ppairs_cnv['px']), set(ppairs_cnv['py'])


# -- Proteasome inhibition data-set
data = read_csv('./data/proteasome_inhibition.csv')[['gene_symbol', 'site_position', 'bortezomib_2hrs_median', 'bortezomib_4hrs_median', 'bortezomib_8hrs_median']]
data = data[[i != '-' for i in data['gene_symbol']]]


# -- Attenuated proteins
cors = read_csv('./tables/proteins_correlations.csv', index_col=0)
att_proteins = set(cors[cors['cluster'] == 1].index)


# -- Plot
plot_df = DataFrame([
    {'gene': data.ix[i, 'gene_symbol'], 'pos': data.ix[i, 'site_position'], 'fc': data.ix[i, c], 'condition': c.split('_')[1]}
        for i in data.index for c in ['bortezomib_2hrs_median', 'bortezomib_4hrs_median', 'bortezomib_8hrs_median']
]).dropna()

plot_df['regulators'] = ['Px' if i in px else ('Py' if i in py else 'Other') for i in plot_df['gene']]
plot_df['attenuation'] = [
    'Not attenuated' if i not in att_proteins else ('Attenuated (> %.1f)' % (np.floor((cors.ix[i, 'diff'] if cors.ix[i, 'diff'] < .5 else .5) * 10) / 10)) for i in plot_df['gene']
]


# Plot regulators
t, pval = ttest_ind(
    plot_df.loc[plot_df['regulators'] == 'Px', 'fc'],
    plot_df.loc[plot_df['regulators'] == 'Py', 'fc'],
    equal_var=False)
print 't: %.2f, p-val: %.2e' % (t, pval)

order, pal, hue_order = ['2hrs', '4hrs', '8hrs'], ['#2980B9', '#E74C3C', '#99A3A4'], ['Px', 'Py', 'Other']

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(plot_df, size=3, aspect=1, legend_out=False)
g = g.map_dataframe(sns.stripplot, x='condition', y='fc', hue='regulators', orient='v', split=True, size=2, jitter=.2, alpha=.6, linewidth=.1, edgecolor='white', order=order, hue_order=hue_order, palette=pal)
g = g.map_dataframe(sns.boxplot, x='condition', y='fc', hue='regulators', orient='v', linewidth=.3, sym='', notch=True, order=order, hue_order=hue_order, palette=pal)
g = g.map(plt.axhline, y=0, ls='-', lw=0.1, c='black', alpha=.5)

g.set_axis_labels('', 'Ubiquitination (Ubq sites log2 FC)')
g.despine(trim=True)
g.add_legend(label_order=hue_order)
plt.suptitle('Proteasome inhibition')
plt.savefig('./reports/proteosome_inhibition_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/proteosome_inhibition_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'


# -- Plot: attenuated proteins
t, pval = ttest_ind(
    plot_df.loc[plot_df['attenuation'] == 'Attenuated', 'fc'],
    plot_df.loc[plot_df['attenuation'] == 'Not attenuated', 'fc'],
    equal_var=False)
print 't: %.2f, p-val: %.2e' % (t, pval)

order = ['2hrs', '4hrs', '8hrs']
hue_order = ['Attenuated (> 0.2)', 'Attenuated (> 0.3)', 'Attenuated (> 0.4)', 'Attenuated (> 0.5)', 'Not attenuated']
pal = dict(zip(*(hue_order, sns.light_palette('#E74C3C', 5).as_hex()[1:] + ['#99A3A4'])))

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(plot_df, size=3.5, aspect=1, legend_out=False)
g = g.map_dataframe(sns.stripplot, x='condition', y='fc', hue='attenuation', orient='v', split=True, size=2, jitter=.2, alpha=.6, linewidth=.1, edgecolor='white', order=order, hue_order=hue_order, palette=pal)
g = g.map_dataframe(sns.boxplot, x='condition', y='fc', hue='attenuation', orient='v', linewidth=.3, sym='', notch=True, order=order, hue_order=hue_order, palette=pal)
g = g.map(plt.axhline, y=0, ls='-', lw=0.1, c='black', alpha=.5)

g.set_axis_labels('', 'Ubiquitination (Ubq sites log2 FC)')
g.despine(trim=True)
g.add_legend(label_order=hue_order)
plt.suptitle('Proteasome inhibition')
plt.savefig('./reports/proteosome_inhibition_attenuation_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/proteosome_inhibition_attenuation_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'
