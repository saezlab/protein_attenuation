import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, AgglomerativeClustering, AffinityPropagation
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from pandas import DataFrame, Series, read_csv
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from cptac import wd, palette, palette_survival, palette_binary


# -- Import
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


# -- Heatmap
plot_df = pancan.dropna().corr(method='pearson')

kmn = AgglomerativeClustering(n_clusters=2, linkage='complete', affinity='correlation').fit(plot_df)
annot_kmn = {k: palette_binary[v] for k, v in dict(zip(*(plot_df.index, kmn.labels_))).items()}

cmap = sns.diverging_palette(220, 10, n=9, as_cmap=True)
sns.set(style='white', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=0.75)
g = sns.clustermap(
    plot_df, metric='correlation', col_colors=[palette[annot.ix[i, 'type']] for i in plot_df],
    row_colors=[annot_kmn[i] for i in plot_df.index], cmap=cmap, figsize=(5, 5),
    xticklabels=False, yticklabels=False
)
plt.suptitle('Pancancer proteomics clustermap')
plt.savefig('%s/reports/pancan_clustermap.png' % wd, bbox_inches='tight', dpi=300)
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

plt.savefig('%s/reports/survival_clustering.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# --
y = Series(dict(zip(*(plot_df.index, kmn.labels_))))
x = pancan.dropna().T.ix[y.index]

lm = LogisticRegression().fit(x, y)
Series(dict(zip(*(x.columns, lm.coef_[0])))).sort_values()
