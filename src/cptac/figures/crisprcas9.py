#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pandas import read_csv, concat


# -- Achilles_v3.3.8
crispr = read_csv('./data/Achilles_v3.3.8.Gs.txt', sep='\t', index_col=0).drop('Description', axis=1)
print crispr.shape


# -- Protein attenuation
p_cor = read_csv('./tables/proteins_correlations.csv', index_col=0)
cor_diff = (p_cor['CNV_Transcriptomics'] - p_cor['CNV_Proteomics']).to_dict()
print p_cor


# -- Regulatory interactions
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')

px, py = set(ppairs_cnv['px']), set(ppairs_cnv['py'])
print len(px), len(py)


# -- Regulators crispr
crispr_df = crispr.unstack().reset_index().copy()
crispr_df.columns = ['cell_line', 'sgrna', 'fc']
crispr_df['tissue'] = [' '.join(i.split('_')[1:]).lower() for i in crispr_df['cell_line']]
crispr_df['type'] = ['px' if i.split('_')[0] in px else ('py' if i.split('_')[0] in py else 'all') for i in crispr_df['sgrna']]
crispr_df['diff'] = [cor_diff[i.split('_')[0]] if i.split('_')[0] in cor_diff else np.nan for i in crispr_df['sgrna']]
print crispr_df

ub, lb = crispr_df['diff'].dropna().quantile(.9), crispr_df['diff'].dropna().quantile(.1)
crispr_df['diff_bin'] = ['attenuated' if i > ub else ('non-attenuated' if i < lb else 'na') for i in crispr_df['diff']]
print crispr_df

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.factorplot(x='fc', y='type', data=crispr_df, kind='box', fliersize=2, linewidth=.3, orient='h', sharex=False, size=1, aspect=1.5, notch=True)
g.map(plt.axvline, x=0, ls='--', lw=0.3, c='black', alpha=.5)
g.despine()
g.set_xlabels('CRISPR-Cas9 (score)')
plt.savefig('./reports/complex_regualtors_crispr_cas9_scores.pdf', bbox_inches='tight')
plt.savefig('./reports/complex_regualtors_crispr_cas9_scores.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.factorplot(x='fc', y='diff_bin', data=crispr_df.dropna(), order=['attenuated', 'na', 'non-attenuated'], kind='box', fliersize=2, linewidth=.3, orient='h', sharex=False, size=1, aspect=1.5, notch=True)
g.map(plt.axvline, x=0, ls='--', lw=0.3, c='black', alpha=.5)
g.despine()
g.set_xlabels('CRISPR-Cas9 (score)')
plt.savefig('./reports/complex_regualtors_crispr_cas9_scores_attenuated.pdf', bbox_inches='tight')
plt.savefig('./reports/complex_regualtors_crispr_cas9_scores_attenuated.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'

