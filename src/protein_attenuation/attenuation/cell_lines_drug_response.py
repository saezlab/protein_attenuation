#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr, ttest_ind
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
from pandas import read_csv, concat, Series, DataFrame, pivot_table
from protein_attenuation.utils import log_likelihood, f_statistic, r_squared


# -- Import data
# Samplesheets
samplesheet = read_csv('./data/sanger_samplesheet.csv', index_col=0).dropna(subset=['TCGA'])
tisues = samplesheet['TCGA'].dropna().reset_index().groupby('TCGA').agg(lambda x: set(x)).to_dict()['sample']

d_samplesheet = read_csv('./data/sanger_drug_samplesheet.csv')
d_targets = d_samplesheet.groupby('Name')['Putative_target'].agg(lambda x: set(x))

# Drug response
drug = read_csv('./data/sanger_drug_response_auc.csv', index_col=1).drop('cosmic', axis=1).T
drug = drug.loc[:, drug.count() > drug.shape[0] * .75]

# Gene expression
trans = read_csv('./data/sanger_gene_experssion_rma.tsv', sep='\t')
trans = pivot_table(trans, index='GENE_NAME', columns='SAMPLE_NAME', values='Z_SCORE', fill_value=np.nan, aggfunc=np.mean)


# -- Import tumour gene-expression signature of protein attenuation
p_attenuation_sig = Series.from_csv('./tables/samples_attenuation_potential_gene_signature.csv')


# -- Cell lines correlation with protein attenuation gene-expression signature
cell_lines_attenuation = trans.ix[p_attenuation_sig.index].corrwith(p_attenuation_sig)
cell_lines_attenuation.name = 'attenuation'


# -- Drug response associations regressions with samples protein attenuation potential
# d = 'AUY922'
def regressions(d):
    df = concat([drug.ix[d], cell_lines_attenuation], axis=1).dropna()

    # Correlation
    cor, pval = pearsonr(df['attenuation'], df[d])

    # Fit model
    lm = LinearRegression().fit(df[['attenuation']], df[d])

    # Predict
    y_true, y_pred = df[d].copy(), Series(dict(zip(*(df.index, lm.predict(df[['attenuation']])))))

    # Log likelihood
    l_lm = log_likelihood(y_true, y_pred)

    # F-statistic
    f, f_pval = f_statistic(y_true, y_pred, len(y_true), df[['attenuation']].shape[1])

    # R-squared
    r = r_squared(y_true, y_pred)

    res = {
        'drug': d, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta': lm.coef_[0], 'pearson': cor, 'pearson_pval': pval
    }
    print '%s: Rsquared: %.2f, F: %.2f, F pval: %.2e, ll: %.2f' % (d, res['rsquared'], res['f'], res['f_pval'], res['ll'])

    return res

ppairs = DataFrame([regressions(d) for d in drug.index])
ppairs['fdr'] = multipletests(ppairs['f_pval'], method='fdr_bh')[1]
ppairs['pearson_fdr'] = multipletests(ppairs['pearson_pval'], method='fdr_bh')[1]
ppairs['targets'] = [';'.join(d_targets.ix[i]) if i in d_targets.index else 'NaN' for i in ppairs['drug']]
ppairs['Putative_target'] = ['; '.join(d_targets[i.split('.')[0]]) if i in d_targets else 'NaN' for i in ppairs['drug']]
ppairs.sort('fdr').to_csv('./tables/drug_response.csv', index=False)
print '[INFO] Drug response associations table: ', './tables/drug_response.csv'


# -- Plot
h_drugs = ['Bortezomib', 'MG-132', 'AUY922', 'SNX-2112', '17-AAG', 'Elesclomol', 'CCT018159', 'Nutlin-3a', 'JNJ-26854165']

# Volcano
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
plt.scatter(
    x=ppairs['beta'], y=-np.log10(ppairs['fdr']), s=25, alpha=.5,
    linewidths=[.5 if i < .05 else 0.01 for i in ppairs['fdr']],
    c=['#E74C3C' if d in h_drugs else ('#99A3A4' if f < .05 else sns.light_palette('#99A3A4').as_hex()[1]) for f, d in ppairs[['fdr', 'drug']].values],
    edgecolor=['black' if f < .05 else sns.light_palette('#99A3A4').as_hex()[1] for f in ppairs['fdr']]
)

# for fdr, beta, d in ppairs[['fdr', 'beta', 'drug']].values:
#     if fdr < .05 and d in h_drugs:
#         plt.text(beta, -np.log10(fdr), '%s' % d, fontsize=6)

plt.axhline(-np.log10(0.01), c='#99A3A4', ls='--', lw=.5, alpha=.7)
plt.axhline(-np.log10(0.05), c='#99A3A4', ls='--', lw=.5, alpha=.7)
plt.axvline(0, c='#99A3A4', ls='-', lw=.3, alpha=.7)

plt.ylim(0)

sns.despine()
plt.ylabel('Adj. p-value (-log10)')
plt.xlabel('Beta')
plt.gcf().set_size_inches(2.5, 5)
plt.savefig('./reports/drug_response_associations_volcano.png', bbox_inches='tight', dpi=600)
plt.savefig('./reports/drug_response_associations_volcano.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Drug response associations volcano: ', './tables/drug_response_associations_volcano.pdf'

# Boxplot
plot_df = ppairs.copy()
plot_df['type'] = ['Proteasome/Chaperone' if i in h_drugs else 'Other' for i in plot_df['drug']]

t, pval = ttest_ind(
    plot_df.loc[plot_df['type'] == 'Other', 'beta'],
    plot_df.loc[plot_df['type'] == 'Proteasome/Chaperone', 'beta'],
    equal_var=False)

hue_order = ['Proteasome/Chaperone', 'Other']
pal = {'Proteasome/Chaperone': '#E74C3C', 'Other': '#99A3A4'}

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(plot_df, size=.75, aspect=3, legend_out=False)
g = g.map_dataframe(sns.stripplot, y='type', x='beta', orient='h', split=True, size=2, jitter=.2, alpha=.6, linewidth=.1, edgecolor='white', order=hue_order, palette=pal)
g = g.map_dataframe(sns.boxplot, y='type', x='beta', orient='h', linewidth=.3, sym='', order=hue_order, palette=pal)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
plt.title('Drug response prediction\n(p-value: %.2e)' % pval)
g.set_axis_labels('Drug response association (beta)', '')
g.despine(trim=True)
plt.savefig('./reports/drug_response_associations_attenuation_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/drug_response_associations_attenuation_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Drug response associations boxplot: ', './tables/drug_response_associations_attenuation_boxplot.pdf'
