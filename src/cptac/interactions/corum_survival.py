import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
import statsmodels.api as sm
from sklearn.mixture.gmm import GMM
from lifelines import KaplanMeierFitter
from matplotlib.gridspec import GridSpec
from lifelines.statistics import logrank_test
from pandas import DataFrame, Series, read_csv
from statsmodels.stats.multitest import multipletests


# -- Imports
# protein complexes scores
pancan = read_csv('%s/tables/pancan_normalised.csv' % wd, index_col=0)
pancan = pancan[pancan.count(1) > (pancan.shape[1] * .5)]
samples = {'-'.join(i.split('-')[:4])[:-1].upper() for i in pancan}
print pancan

# correlating pairs
p_pairs = read_csv('%s/tables/top_correlated_protein_pairs.csv' % wd)
p_pairs = p_pairs[p_pairs['score'] > .5]
print p_pairs

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0).ix[samples, ['DAYS_TO_LAST_FOLLOWUP', 'VITAL_STATUS']]
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
print clinical


# --
# p1, p2 = 'ACOT2', 'ACOT1'
def gmm_survival(p1, p2):
    name = '%s-%s' % (p1, p2)
    print name

    # -- Residuals
    # Get protein measurements
    x = pancan.ix[[p1, p2]].T.dropna()[p1]
    y = pancan.ix[[p1, p2]].T.dropna()[p2]

    lm = sm.OLS(y, sm.add_constant(x)).fit()
    residuals = DataFrame(lm.resid)

    # -- GMM
    gmm = GMM(n_components=3).fit(residuals)

    groups = {
        'low': gmm.means_.argmin(),
        'high': gmm.means_.argmax(),
        'neutral': list({0, 1, 2} - {gmm.means_.argmin(), gmm.means_.argmax()})[0]
    }

    groups = {g: set(residuals[gmm.predict(residuals) == groups[g]].index) for g in groups}
    groups_ = {g: {'-'.join(i.split('-')[:4])[:-1].upper() for i in groups[g]} for g in groups}

    logrank = np.nan

    if len(groups['low']) > 3 and len(groups['high']) > 3 and len(groups['neutral']) > 3 and len(residuals) > (.5 * pancan.shape[1]):
        # -- Logrank test
        logrank = logrank_test(
            clinical.ix[groups_['high'], 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[groups_['low'], 'DAYS_TO_LAST_FOLLOWUP'],
            clinical.ix[groups_['high'], 'VITAL_STATUS'], clinical.ix[groups_['low'], 'VITAL_STATUS']
        )

        # -- Plot
        if logrank.p_value < .001:
            sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
            fig, gs = plt.figure(figsize=(10, 3)), GridSpec(1, 3, hspace=.3)
            palette = {'high': '#2ecc71', 'neutral': '#808080', 'low': '#e74c3c'}

            # Survival curves
            ax = plt.subplot(gs[0])

            kmf = KaplanMeierFitter()
            for g in groups:
                kmf.fit(clinical.ix[groups_[g], 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[groups_[g], 'VITAL_STATUS'], label=g)
                kmf.plot(ax=ax, ci_force_lines=False, color=palette[g])

            ax.set_xlabel('Timeline (days)')
            ax.set_ylabel('Survival probability')
            sns.despine(ax=ax)

            # Histogram
            ax = plt.subplot(gs[1])

            for g in groups:
                sns.distplot(residuals.ix[groups[g]], color=palette[g], kde=False, ax=ax)

            sns.despine(ax=ax)
            ax.set_xlabel('Protein-pair residuals')
            ax.set_ylabel('Counts')

            # Scatter
            ax = plt.subplot(gs[2])

            for g in groups:
                sns.regplot(x.ix[groups[g]], y.ix[groups[g]], ax=ax, label=g, color=palette[g], fit_reg=False)

            ax.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)
            ax.axvline(0, ls='--', lw=0.1, c='black', alpha=.3)
            sns.despine(ax=ax)

            # Export
            # fig.suptitle('%s; p-value: %.2e' % (name, logrank.p_value))
            plt.savefig('%s/reports/survival_complex_%s.pdf' % (wd, name), bbox_inches='tight')
            plt.close('all')

            print 'interaction: %s; p-value: %.2e, samples: %d' % (name, logrank.p_value, len(residuals))

        return {'pval': logrank.p_value}

p_survival = DataFrame({'%s - %s' % (p1, p2): gmm_survival(p1, p2) for p1, p2 in p_pairs[['p1', 'p2']].values}).T.dropna()
p_survival['fdr'] = multipletests(p_survival['pval'], method='fdr_bh')[1]
print p_survival.sort('fdr')

