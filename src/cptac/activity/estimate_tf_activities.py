#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, default_color, palette_cnv_number
from cptac.utils import ztest
from pandas import DataFrame, read_csv
from statsmodels.stats.multitest import multipletests


# -- TFs regulons
regulons = read_csv('%s/files/tfs_regulons.csv' % wd)
regulons = regulons.groupby('set')['gene'].agg(lambda x: set(x)).to_dict()
regulons = {k: regulons[k].difference(set(k)) for k in regulons if len(regulons[k].difference(set(k))) > 0}
print len(regulons)


# -- Transcriptomics
trans = read_csv('%s/data/tcga_rnaseq_fpkm.csv' % wd, index_col=0).replace(0, np.nan)
trans = np.log2(trans)
print trans


# -- Estimate TF activity
# s, t = 'TCGA-E2-A15A-01', 'ZEB1'
# targets, mu, var = trans.ix[regulons[t], s].dropna(), trans[s].mean(), trans[s].var()
tf_activity = [[s, t] + list(ztest(trans.ix[regulons[t], s].dropna(), trans[s].mean(), trans[s].var())) for s in trans for t in regulons if len(trans.ix[regulons[t], s].dropna()) > 0]
tf_activity = DataFrame(tf_activity, columns=['sample', 'tf', 'z', 'pval', 'mean', 'targets'])
tf_activity['FDR'] = multipletests(tf_activity['pval'],  method='fdr_bh')[1]
tf_activity.to_csv('%s/tables/tf_activities.csv' % wd)
print tf_activity.sort('pval')


# # -- QQ-plot
# sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
#
# plot_df = DataFrame({
#     'x': sorted([-np.log10(np.float(i) / len(tf_activity)) for i in np.arange(1, len(tf_activity) + 1)]),
#     'y': sorted(-np.log10(tf_activity['pval']))
# })
#
# g = sns.regplot('x', 'y', plot_df, fit_reg=False, ci=False, color=default_color, line_kws={'lw': .3})
# g.set_xlim(0)
# g.set_ylim(0)
# plt.plot(plt.xlim(), plt.xlim(), 'k--', lw=.3)
# plt.xlabel('Theoretical -log(P)')
# plt.ylabel('Observed -log(P)')
# plt.title('TF activities')
# sns.despine(trim=True)
# plt.gcf().set_size_inches(3, 3)
# plt.savefig('%s/reports/tf_activities_qqplot.pdf' % wd, bbox_inches='tight')
# plt.close('all')
# print '[INFO] Plot done'


# -- CNV validation
# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
tf_activity['cnv'] = [cnv.ix[t, s] if t in cnv.index and s in cnv.columns else np.nan for s, t in tf_activity[['sample', 'tf']].values]

sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
plt.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)
sns.violinplot('cnv', 'z', data=tf_activity, palette=palette_cnv_number, linewidth=.3, cut=0, inner='quartiles')
sns.despine(trim=True)
plt.gcf().set_size_inches(3, 3)
plt.savefig('%s/reports/tf_activities_boxplots.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
