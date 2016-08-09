import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from cptac import wd, palette, palette_cnv_number
from matplotlib.gridspec import GridSpec
from sklearn.linear_model import LinearRegression
from matplotlib_venn import venn3, venn3_circles
from statsmodels.stats.multitest import multipletests
from lifelines import KaplanMeierFitter
from sklearn.metrics.ranking import roc_curve, auc
from lifelines.statistics import logrank_test
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pandas import DataFrame, Series, read_csv, concat, merge, pivot_table
from pymist.utils.map_peptide_sequence import read_uniprot_genename, read_fasta


# -- Imports
# Protein complexes activities
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print proteomics

# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print 'transcriptomics', transcriptomics.shape

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0).dropna(subset=['VITAL_STATUS', 'DAYS_TO_LAST_FOLLOWUP'])
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
clinical = clinical[clinical['DAYS_TO_LAST_FOLLOWUP'] > 1]
print clinical

# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
print 'cnv', cnv.shape

# Complexes activity
p_complex = read_csv('%s/tables/protein_complexes_activities.tsv' % wd, index_col=0, sep='\t')
p_complex = pivot_table(p_complex, index='complex', columns='sample', values='z')
print p_complex

# P-pairs
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot} for k, v in get_complexes_dict().items()}
corum_n = get_complexes_name()

ppairs = read_csv('%s/tables/ppairs_cnv_regulation.csv' % wd)
ppairs_corum = {i: [(c, corum_n[c]) for c in corum if i in corum[c]] for i in ppairs['px']}
ppairs_signif = set(ppairs['px']).union(ppairs['py'])
ppairs_complexes = {c for px, py in ppairs[['px', 'py']].values for c in corum if px in corum[c] and py in corum[c]}
print ppairs


# -- Overlap
samples = set(proteomics).intersection(clinical.index).intersection(p_complex).intersection(transcriptomics)
print len(samples)


# p = 'AP1M1'
def protein_residual(p):
    y = proteomics.ix[p, samples].dropna()
    x = transcriptomics.ix[[p], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[p] - lm.intercept_

    return y_
residuals = DataFrame({p: protein_residual(p) for p in set(ppairs['px']).union(ppairs['py'])}).T
print 'residuals', residuals.shape


# -- Survival
# c = 351
def survival(c):
    y = p_complex.ix[c, samples].dropna()

    samples_up = set(y[y > np.percentile(y, 85)].index)
    samples_down = set(y[y < np.percentile(y, 15)].index)

    print corum_n[c], len(samples_down), len(samples_up)

    if len(samples_down) > 0 and len(samples_up) > 0:
        logrank = logrank_test(
            clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'],
            clinical.ix[samples_up, 'VITAL_STATUS'], clinical.ix[samples_down, 'VITAL_STATUS']
        )

        res = {
            'complex': c,
            'ttest': logrank.test_statistic, 'logrank': logrank.p_value,
        }
        print res

        return res

res = [survival(c) for c in ppairs_complexes]
res = DataFrame([i for i in res if i])
res['fdr'] = multipletests(res['logrank'], method='fdr_bh')[1]
print res.sort('fdr')


# --
plot_df = [938, 936]

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs, pos = GridSpec(len(plot_df), 1, hspace=.5), 0
for c in plot_df:
    # Define survival test
    y = p_complex.ix[c, samples].dropna()
    assoc_pairs = ['Py (%s) ~ Px (%s)' % (py, px) for px, py in ppairs[['px', 'py']].values if px in corum[c] and py in corum[c]]

    samples_up = set(y[y > np.percentile(y, 85)].index)
    samples_down = set(y[y < np.percentile(y, 15)].index)
    samples_bkg = samples.difference(samples_up.union(samples_down))

    logrank = logrank_test(
        clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'],
        clinical.ix[samples_up, 'VITAL_STATUS'], clinical.ix[samples_down, 'VITAL_STATUS']
    )

    # Plot
    ax = plt.subplot(gs[pos])

    kmf = KaplanMeierFitter()

    kmf.fit(clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_up, 'VITAL_STATUS'], label='Up')
    kmf.plot(ax=ax, ci_force_lines=False, color=palette_cnv_number[1])

    kmf.fit(clinical.ix[samples_bkg, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_bkg, 'VITAL_STATUS'], label='Background')
    kmf.plot(ax=ax, ci_force_lines=False, color=palette_cnv_number[0])

    kmf.fit(clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_down, 'VITAL_STATUS'], label='Down')
    kmf.plot(ax=ax, ci_force_lines=False, color=palette_cnv_number[-1])

    sns.despine(ax=ax)

    ax.set_title('%s\nLogrank test: %.2e\n%s' % (corum_n[c], logrank.p_value, '; '.join(assoc_pairs)))
    ax.set_xlabel('Timeline (days)')
    ax.set_ylabel('Survival fraction')

    pos += 1

plt.gcf().set_size_inches(4, 3 * len(plot_df))
plt.savefig('%s/reports/survival_protein_complexes.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# -- Survival - CNV
# p = 'CAPZB'
def survival(p):
    y = residuals.ix[p, samples].dropna()

    # samples_up = set(y[y > 0].index)
    # samples_down = set(y[y < 0].index)
    samples_up = set(y[y > np.percentile(y, 85)].index)
    samples_down = set(y[y < np.percentile(y, 15)].index)

    print p, len(samples_down), len(samples_up)

    if len(samples_down) > 0 and len(samples_up) > 0:
        logrank = logrank_test(
            clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'],
            clinical.ix[samples_up, 'VITAL_STATUS'], clinical.ix[samples_down, 'VITAL_STATUS']
        )

        res = {
            'px': p,
            'ttest': logrank.test_statistic, 'logrank': logrank.p_value,
        }
        print res

        return res

res_cnv = [survival(c) for c in set(ppairs['py'])]
res_cnv = DataFrame([i for i in res_cnv if i])
res_cnv['fdr'] = multipletests(res_cnv['logrank'], method='fdr_bh')[1]
print res_cnv.sort('fdr')


# --
plot_df = ['SSRP1']

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs, pos = GridSpec(len(plot_df), 1, hspace=.5), 0
for p in plot_df:
    y = residuals.ix[p, samples].dropna()

    samples_up = set(y[y > np.percentile(y, 85)].index)
    samples_down = set(y[y < np.percentile(y, 15)].index)
    samples_bkg = samples.difference(samples_up.union(samples_down))

    logrank = logrank_test(
        clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'],
        clinical.ix[samples_up, 'VITAL_STATUS'], clinical.ix[samples_down, 'VITAL_STATUS']
    )

    # Plot
    ax = plt.subplot(gs[pos])

    kmf = KaplanMeierFitter()

    kmf.fit(clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_up, 'VITAL_STATUS'], label='Up')
    kmf.plot(ax=ax, ci_force_lines=False, color=palette_cnv_number[1])

    kmf.fit(clinical.ix[samples_bkg, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_bkg, 'VITAL_STATUS'], label='Background')
    kmf.plot(ax=ax, ci_force_lines=False, color=palette_cnv_number[0])

    kmf.fit(clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_down, 'VITAL_STATUS'], label='Down')
    kmf.plot(ax=ax, ci_force_lines=False, color=palette_cnv_number[-1])

    sns.despine(ax=ax)

    ax.set_title('%s\nLogrank test: %.2e' % (p, logrank.p_value))
    ax.set_xlabel('Timeline (days)')
    ax.set_ylabel('Survival fraction')

    pos += 1

plt.gcf().set_size_inches(4, 3 * len(plot_df))
plt.savefig('%s/reports/survival_protein_cnv.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
