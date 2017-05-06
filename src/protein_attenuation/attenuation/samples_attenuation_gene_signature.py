#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from protein_attenuation.enrichment.gsea import gsea
from sklearn.feature_selection import SelectFdr, f_classif
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.linear_model.logistic import LogisticRegressionCV
from statsmodels.stats.multitest import multipletests
from protein_attenuation.utils import read_gmt
from pandas import DataFrame, Series, read_csv, concat


# -- Imports
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
s_attenuation = read_csv('./tables/sample_attenuation_table.csv', index_col=0)
signature = Series.from_csv('./tables/samples_attenuation_potential_gene_signature.csv')

# -- GO terms
msigdb_go_bp = read_gmt('./files/c5.bp.v5.1.symbols.gmt')
msigdb_go_cc = read_gmt('./files/c5.cc.v5.1.symbols.gmt')
msigdb_go_mf = read_gmt('./files/c5.mf.v5.1.symbols.gmt')


# -- Logistic regression with bootstrap
x, y = transcriptomics[s_attenuation.index].T, np.logical_not(s_attenuation['cluster']).astype(int)

scores, cv = [], StratifiedShuffleSplit(n_splits=1000, test_size=.3)
for train, test in cv.split(x, y):
    xs, ys = x.ix[train].copy(), y.ix[train].copy()

    fs = SelectFdr(f_classif).fit(xs, ys)
    b_features = list(xs.loc[:, fs.pvalues_ < 0.05].columns)

    lm = LogisticRegressionCV().fit(xs[b_features], ys)

    scores.append(roc_auc_score(y.ix[test], lm.decision_function(x.ix[test, b_features])))

Series(scores).to_csv('./tables/samples_attenuation_logistic_aucs.csv')
print(np.mean(scores), np.percentile(scores, 95), np.percentile(scores, 5))

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'lines.linewidth': .75})
sns.distplot(scores, color='#34495e')
for v in [np.mean(scores), np.percentile(scores, 95), np.percentile(scores, 5)]:
    plt.axvline(v, ls='--', lw=.3, color='#e74c3c')
plt.axvline(.5, ls='-', lw=.3, color='#dfdfdf')
sns.despine()
plt.title('Logistic classification (mean AROC: %.2f)' % np.mean(scores))
plt.xlabel('AROC')
plt.ylabel('Density')
plt.gcf().set_size_inches(5, 3)
plt.savefig('./reports/samples_attenuation_logistic_aucs_histogram.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Gene signature enrichment analysis
dataset = signature.to_dict()
signatures = {'BP': msigdb_go_bp, 'CC': msigdb_go_cc, 'MF': msigdb_go_mf}

df_enrichment = [(t, sig, len(db[sig].intersection(dataset)), gsea(dataset, db[sig], 1000)) for t, db in signatures.items() for sig in db]
df_enrichment = DataFrame([{'type': t, 'signature': s, 'length': l, 'escore': es, 'pvalue': pval} for t, s, l, (es, pval) in df_enrichment]).dropna()
df_enrichment.sort(['pvalue', 'escore']).to_csv('./tables/enrichment_sample_attenuation_gene_signature.csv', index=False)
# df_enrichment = read_csv('./tables/enrichment_sample_attenuation_gene_signature.csv')
print(df_enrichment)

# Plot - Striplot
plot_df = df_enrichment[(df_enrichment['length'] > 5)].copy()
plot_df['fdr'] = multipletests(plot_df['pvalue'], method='fdr_bh')[1]
plot_df = plot_df[plot_df['fdr'].abs() < .05].sort('escore')
plot_df['signature'] = [i.replace('_', ' ').lower() for i in plot_df['signature']]
plot_df = concat([plot_df.head(30), plot_df.tail(30)])

pal = dict(zip(*(set(plot_df['type']), sns.color_palette('Set1', n_colors=4).as_hex())))

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
sns.stripplot(y='signature', x='escore', hue='type', data=plot_df, palette=pal, s=8)
plt.axvline(0, ls='-', lw=0.3, c='black', alpha=.5)
sns.despine(trim=True)
plt.xlabel('Enrichment score (GSEA)')
plt.ylabel('')
plt.title('Attenuation potential - gene signature')
plt.gcf().set_size_inches(2, 10)
plt.savefig('./reports/sample_attenuation_enrichment.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/sample_attenuation_enrichment.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Protein attenuation complexes: ', './reports/protein_correlation_difference_enrichment.pdf'
