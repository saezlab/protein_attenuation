import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from sklearn.mixture.gmm import GMM
from lifelines import KaplanMeierFitter
from matplotlib.gridspec import GridSpec
from lifelines.statistics import logrank_test
from pandas import DataFrame, Series, read_csv
from statsmodels.stats.multitest import multipletests


# -- Imports
# protein complexes scores
corum = read_csv('%s/tables/corum_scores.csv' % wd, index_col=0)
samples = {'-'.join(i.split('-')[:4])[:-1].upper() for i in corum}
corum = corum[corum.count(1) > (corum.shape[1] * .5)]

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0).ix[samples, ['DAYS_TO_LAST_FOLLOWUP', 'VITAL_STATUS']]
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]


# --
# c = 'APP-TIMM23 complex'
def complex_gmm(c):
    # -- Get complex scores
    complex_scores = corum.ix[[c]].T.dropna()

    # -- GMM
    gmm = GMM(n_components=2).fit(complex_scores)

    group0 = set(complex_scores[np.logical_not(gmm.predict(complex_scores).astype(np.bool))].index)
    group0_ = {'-'.join(i.split('-')[:4])[:-1].upper() for i in group0}

    group1 = set(complex_scores[gmm.predict(complex_scores).astype(np.bool)].index)
    group1_ = {'-'.join(i.split('-')[:4])[:-1].upper() for i in group1}

    # -- Logrank test
    logrank = logrank_test(
        clinical.ix[group0_, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[group1_, 'DAYS_TO_LAST_FOLLOWUP'],
        clinical.ix[group0_, 'VITAL_STATUS'], clinical.ix[group1_, 'VITAL_STATUS']
    )

    # -- Plot
    if logrank.p_value < .001:
        sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
        fig, gs = plt.figure(figsize=(7, 3)), GridSpec(1, 2, hspace=.3)
        palette = {'Group 0': '#808080', 'Group 1': '#e74c3c'}

        # Survival curves
        ax = plt.subplot(gs[0])

        kmf = KaplanMeierFitter()
        for n, g in [('Group 0', group0_), ('Group 1', group1_)]:
            kmf.fit(clinical.ix[g, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[g, 'VITAL_STATUS'], label=n)
            kmf.plot(ax=ax, ci_force_lines=False, color=palette[n])

        ax.set_xlabel('Timeline (days)')
        ax.set_ylabel('Survival probability')
        sns.despine(ax=ax)

        # Histogram
        ax = plt.subplot(gs[1])

        sns.distplot(complex_scores.ix[group1], color=palette['Group 1'], kde=False, ax=ax)
        sns.distplot(complex_scores.ix[group0], color=palette['Group 0'], kde=False, ax=ax)

        sns.despine(ax=ax)
        ax.set_xlabel('Protein complex expression score')
        ax.set_ylabel('Counts')

        # Export
        fig.suptitle('Complex: %s; p-value: %.2e' % (c, logrank.p_value))
        plt.savefig('%s/reports/survival_complex_%s.pdf' % (wd, c.replace('/', '_')), bbox_inches='tight')
        plt.close('all')

        print 'complex: %s; p-value: %.2e, samples: %d' % (c, logrank.p_value, len(complex_scores))

    return {'complex': c, 'pval': logrank.p_value}

corum_gmm = DataFrame({c: complex_gmm(c) for c in corum.index}).T
corum_gmm['fdr'] = multipletests(corum_gmm['pval'], method='fdr_bh')[1]
print corum_gmm.sort('fdr')

