import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from cptac import wd, palette, palette_binary
from pandas import DataFrame, Series, read_csv, pivot_table
from pandas.stats.misc import zscore
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from sklearn.linear_model import LogisticRegression
from sklearn.cluster.hierarchical import AgglomerativeClustering
from pymist.utils.corumdb import get_complexes_dict, get_complexes_name
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# proteomics
pancan = read_csv('%s/tables/pancan_preprocessed_normalised.csv' % wd, index_col=0)
pancan = pancan[pancan.count(1) > (pancan.shape[1] * .5)]
annot = read_csv('%s/tables/samplesheet.csv' % wd, index_col=0)
print pancan

# clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0)
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
clinical = clinical[clinical['DAYS_TO_LAST_FOLLOWUP'] > 1]
print clinical

# cnv
cnv_df = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
cnv_df['barcode'] = [i[:15] for i in cnv_df['barcode']]
cnv = pivot_table(cnv_df, index=['hgnc', 'cytoband'], columns='barcode', values='gistic', fill_value=0).astype(np.int)
print cnv

# rnaseq
gexp = read_csv('%s/data/tcga_rnaseq.tsv' % wd, sep='\t', index_col=0)
gexp = zscore(gexp.T).T
print gexp

# corum
uniprot = read_uniprot_genename()
corum_dict = {k: {uniprot[p][0]: 1 for p in v if p in uniprot and uniprot[p][0] in pancan.index} for k, v in get_complexes_dict().items()}
corum_dict = {k: corum_dict[k] for k in corum_dict if len(corum_dict[k]) > 1}
corum = DataFrame(corum_dict).replace(np.nan, 0).astype(np.int)
corum_names = get_complexes_name()
print corum

# linear regression
# s = 'TCGA-E2-A15A-01A-41-A21W-30'
c_scores = {}
for s in pancan:
    y = pancan.ix[corum.index, s].dropna()

    x = corum.ix[y.index]
    x = x.loc[:, x.sum() > 1]

    lm = sm.OLS(y, x).fit()
    print s, lm.rsquared

    c_scores[s] = lm.tvalues

c_scores = DataFrame(c_scores)
c_scores.to_csv('%s/tables/complexes_scores.csv' % wd)
print c_scores

#
plot_df = c_scores.loc[:, c_scores.std() > 1].dropna().corr(method='pearson')

kmn = AgglomerativeClustering(n_clusters=2, linkage='complete', affinity='correlation').fit(plot_df)
annot_kmn = {k: palette_binary[v] for k, v in dict(zip(*(plot_df.index, kmn.labels_))).items()}

cmap = sns.diverging_palette(220, 10, n=9, as_cmap=True)
sns.set(style='white', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=0.75)
g = sns.clustermap(
    plot_df, metric='correlation', col_colors=[palette[annot.ix[i, 'type']] for i in plot_df],
    row_colors=[annot_kmn[i] for i in plot_df.index], cmap=cmap, figsize=(5, 5),
    xticklabels=False, yticklabels=False
)
plt.suptitle('Pancancer complexes clustermap')
plt.savefig('%s/reports/pancan_clustermap_complexes.png' % wd, bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Plot done'


# -- Survival
groups = DataFrame(Series(dict(zip(*(plot_df.index, kmn.labels_)))), columns=['cluster'])
groups = groups[[i[:15] in clinical.index for i in groups.index]]
groups['days'] = [clinical.ix[i[:15], 'DAYS_TO_LAST_FOLLOWUP'] for i in groups.index]
groups['status'] = [clinical.ix[i[:15], 'VITAL_STATUS'] for i in groups.index]

logrank = logrank_test(
    groups.loc[groups['cluster'] == 0, 'days'], groups.loc[groups['cluster'] == 1, 'days'],
    groups.loc[groups['cluster'] == 0, 'status'], groups.loc[groups['cluster'] == 1, 'status']
)
print logrank

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
ax = plt.subplot(111)

kmf = KaplanMeierFitter()
for g in [0, 1]:
    kmf.fit(
        groups.loc[groups['cluster'] == g, 'days'],
        event_observed=groups.loc[groups['cluster'] == g, 'status'],
        label=str(g)
    )
    kmf.plot(ci_force_lines=False, color=palette_binary[g], ax=ax)

plt.xlabel('Timeline (days)')
plt.ylabel('Survival probability')
sns.despine()

plt.savefig('%s/reports/survival_clustering_complexes.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# --
y = Series(dict(zip(*(plot_df.index, kmn.labels_))))
x = c_scores[y.index].dropna().T

lm = LogisticRegression().fit(x, y)
coefs = DataFrame(Series(dict(zip(*(x.columns, lm.coef_[0])))), columns=['coef'])
coefs['name'] = [corum_names[i] for i in coefs.index]
print coefs.sort('coef')


# --
cnv_down = {(b, g) for b, g in cnv_df[cnv_df['gistic'] <= -2][['barcode', 'hgnc']].values}
cnv_up = {(b, g) for b, g in cnv_df[cnv_df['gistic'] >= 2][['barcode', 'hgnc']].values}


plot_df = pancan.unstack().reset_index().dropna()
plot_df.columns = ['barcode', 'protein', 'score']
plot_df['cnv'] = [-1 if (b[:15], p) in cnv_down else (1 if (b[:15], p) in cnv_up else 0) for b, p in plot_df[['barcode', 'protein']].values]
print plot_df

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.boxplot(x='cnv', y='score', data=plot_df)
sns.despine(bottom=True)
plt.gcf().set_size_inches(3, 5)
plt.savefig('%s/reports/pancan_protein_cnv_boxplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


plot_df = gexp.unstack().reset_index().dropna()
plot_df.columns = ['barcode', 'protein', 'score']
plot_df['cnv'] = [-1 if (b[:15], p) in cnv_down else (1 if (b[:15], p) in cnv_up else 0) for b, p in plot_df[['barcode', 'protein']].values]
print plot_df

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.boxplot(x='cnv', y='score', data=plot_df)
sns.despine(bottom=True)
plt.gcf().set_size_inches(3, 5)
plt.savefig('%s/reports/pancan_gexp_cnv_boxplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
