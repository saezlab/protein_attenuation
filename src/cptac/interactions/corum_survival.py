import igraph
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import matplotlib.pyplot as plt
from sklearn.mixture.gmm import GMM
from cptac import wd
from lifelines import KaplanMeierFitter
from matplotlib.gridspec import GridSpec
from lifelines.statistics import logrank_test
from pandas import DataFrame, Series, read_csv
from pymist.utils.corumdb import get_complexes_dict
from statsmodels.stats.multitest import multipletests
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
# Uniprot
uniprot = read_uniprot_genename()

# CORUM
corum = get_complexes_dict()
corum = {k: {uniprot[i][0] for i in v if i in uniprot} for k, v in corum.items()}

# protein complexes scores
c_activity = read_csv('%s/tables/complex_activity.csv' % wd, index_col=0)
c_activity = c_activity[c_activity.count(1) > (c_activity.shape[1] * .25)]
samples = {'-'.join(i.split('-')[:4])[:-1].upper() for i in c_activity}
print c_activity

# correlating pairs
p_pairs = read_csv('%s/tables/top_correlated_protein_pairs.csv' % wd)
p_pairs = p_pairs[p_pairs['corum'] == 1]
print p_pairs

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0).ix[samples, ['DAYS_TO_LAST_FOLLOWUP', 'VITAL_STATUS']]
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
print clinical


# --
# p1, p2 = 'ACOT2', 'ACOT1'
# c = 'SNARE complex (VAMP2, SNAP25, STX1a, STX3, CPLX1, CPLX3, CPLX4)'
def gmm_survival(c):
    print c

    # -- Residuals
    # Get protein measurements
    scores = DataFrame(c_activity.ix[c].dropna().abs())

    # -- GMM
    gmm = GMM(n_components=2).fit(scores)

    groups = {
        'low': gmm.means_.argmin(),
        'high': gmm.means_.argmax()
    }

    groups = {g: set(scores[gmm.predict(scores) == groups[g]].index) for g in groups}
    groups_ = {g: {'-'.join(i.split('-')[:4])[:-1].upper() for i in groups[g]} for g in groups}

    if len(groups['low']) > 1 and len(groups['high']) > 1:
        # -- Logrank test
        logrank = logrank_test(
            clinical.ix[groups_['high'], 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[groups_['low'], 'DAYS_TO_LAST_FOLLOWUP'],
            clinical.ix[groups_['high'], 'VITAL_STATUS'], clinical.ix[groups_['low'], 'VITAL_STATUS']
        )

        # -- Plot
        if logrank.p_value < .001:
            sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
            fig, gs = plt.figure(figsize=(7, 3)), GridSpec(1, 2, hspace=.3)
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
                sns.distplot(scores.ix[groups[g]], color=palette[g], kde=False, ax=ax)

            sns.despine(ax=ax)
            ax.set_xlabel('Protein-pair residuals')
            ax.set_ylabel('Counts')

            # # Scatter
            # ax = plt.subplot(gs[2])
            #
            # for g in groups:
            #     sns.regplot(x.ix[groups[g]], y.ix[groups[g]], ax=ax, label=g, color=palette[g], fit_reg=False)
            #
            # ax.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)
            # ax.axvline(0, ls='--', lw=0.1, c='black', alpha=.3)
            # sns.despine(ax=ax)

            # Export
            # fig.suptitle('%s; p-value: %.2e' % (name, logrank.p_value))
            plt.savefig('%s/reports/survival_complex_%s.pdf' % (wd, c.replace('/', '_')), bbox_inches='tight')
            plt.close('all')

            print 'interaction: %s; p-value: %.2e, samples: %d' % (c, logrank.p_value, len(c))

        return {'pval': logrank.p_value}

    return {'pval': np.nan}

p_survival = DataFrame({c: gmm_survival(c) for c in c_activity.index}).T.dropna()
p_survival['fdr'] = multipletests(p_survival['pval'], method='fdr_bh')[1]
print p_survival.sort('fdr')

