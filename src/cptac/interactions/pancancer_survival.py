import igraph
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import matplotlib.pyplot as plt
from sklearn.mixture.gmm import GMM
from cptac import wd, palette_survival
from lifelines import KaplanMeierFitter
from matplotlib.gridspec import GridSpec
from lifelines.statistics import logrank_test
from pandas import DataFrame, Series, read_csv
from pymist.utils.corumdb import get_complexes_dict
from statsmodels.stats.multitest import multipletests
from sklearn.covariance.outlier_detection import EllipticEnvelope
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
# clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0)
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
clinical = clinical[clinical['DAYS_TO_LAST_FOLLOWUP'] > 1]
print clinical

# pancancer proteomics
pancan = read_csv('%s/tables/pancan_preprocessed_normalised.csv' % wd, index_col=0)
pancan = pancan[pancan.count(1) > (pancan.shape[1] * .5)]
pancan = pancan.loc[:, ['-'.join(i.split('-')[:4])[:-1].upper() in clinical.index for i in pancan]]
print pancan

# Top correlated proteins
pps = read_csv('%s/tables/top_correlated_protein_pairs.csv' % wd)
print pps


# --
# p1, p2 = 'ACTR3', 'ACTR2'
# c = 'SNARE complex (VAMP2, SNAP25, STX1a, STX3, CPLX1, CPLX3, CPLX4)'
# c = 'ABCB8'
def gmm_survival(p1, p2):
    name = '%s-%s' % (p1, p2)
    print name

    # -- Residuals
    # Get protein measurements
    x = pancan.ix[[p1, p2]].T.dropna()[p1]
    y = pancan.ix[[p1, p2]].T.dropna()[p2]

    lm = sm.OLS(y, sm.add_constant(x)).fit()

    # -- Outliers
    groups = {'low': set(lm.resid.ix[lm.resid.abs() <= 1].index), 'high': set(lm.resid.ix[lm.resid.abs() > 1.].index)}
    groups_ = {g: {'-'.join(i.split('-')[:4])[:-1].upper() for i in groups[g]} for g in groups}

    if len(groups['low']) > 1 and len(groups['high']) > 1:
        # -- Logrank test
        logrank = logrank_test(
            clinical.ix[groups_['high'], 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[groups_['low'], 'DAYS_TO_LAST_FOLLOWUP'],
            clinical.ix[groups_['high'], 'VITAL_STATUS'], clinical.ix[groups_['low'], 'VITAL_STATUS']
        )

        # -- Plot
        if logrank.p_value < .001 and len(groups['low']) >= 5 and len(groups['high']) >= 5:
            sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
            fig, gs = plt.figure(figsize=(10, 3)), GridSpec(1, 3, hspace=.3)

            # Survival curves
            ax = plt.subplot(gs[0])

            kmf = KaplanMeierFitter()
            for g in groups:
                kmf.fit(clinical.ix[groups_[g], 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[groups_[g], 'VITAL_STATUS'], label=g)
                kmf.plot(ax=ax, ci_force_lines=False, color=palette_survival[g])

            ax.set_xlabel('Timeline (days)')
            ax.set_ylabel('Survival probability')
            sns.despine(ax=ax)

            # Histogram
            ax = plt.subplot(gs[1])

            for g in groups:
                sns.distplot(lm.resid.ix[groups[g]], color=palette_survival[g], kde=False, ax=ax)

            sns.despine(ax=ax)
            ax.set_xlabel('Protein-pair residuals')
            ax.set_ylabel('Counts')

            # Scatter
            ax = plt.subplot(gs[2])

            for g in groups:
                sns.regplot(x.ix[groups[g]], y.ix[groups[g]], ax=ax, label=g, color=palette_survival[g], fit_reg=False)

            ax.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)
            ax.axvline(0, ls='--', lw=0.1, c='black', alpha=.3)
            sns.despine(ax=ax)

            # Export
            # fig.suptitle('%s; p-value: %.2e' % (name, logrank.p_value))
            plt.savefig('%s/reports/survival_complex_%s.pdf' % (wd, name), bbox_inches='tight')
            plt.close('all')
            print 'interaction: %s; p-value: %.2e, samples: %d' % (name, logrank.p_value, len(lm.resid))

        return {'pval': logrank.p_value, 'low': len(groups['low']), 'high': len(groups['high'])}

    return {'pval': np.nan}

p_survival = DataFrame({'%s - %s' % (p1, p2): gmm_survival(p1, p2) for p1, p2 in pps[['p1', 'p2']].values}).T.dropna()
p_survival['fdr'] = multipletests(p_survival['pval'], method='fdr_bh')[1]
print p_survival.sort('fdr')

