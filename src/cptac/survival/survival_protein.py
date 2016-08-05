import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from cptac import wd, palette, palette_cnv_number
from matplotlib.gridspec import GridSpec
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

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0).dropna(subset=['VITAL_STATUS', 'DAYS_TO_LAST_FOLLOWUP'])
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
clinical = clinical[clinical['DAYS_TO_LAST_FOLLOWUP'] > 1]
print clinical

# Complexes activity
p_complex = read_csv('%s/tables/protein_complexes_activities_mean.tsv' % wd, index_col=0, sep='\t')
print p_complex

# P-pairs
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot} for k, v in get_complexes_dict().items()}
corum_n = get_complexes_name()

ppairs = read_csv('%s/tables/ppairs_cnv_regulation.csv' % wd)
ppairs_corum = {i: [(c, corum_n[c]) for c in corum if i in corum[c]] for i in ppairs['px']}
ppairs_signif = set(ppairs['px']).union(ppairs['py'])
print ppairs


# -- Overlap
samples = set(proteomics).intersection(clinical.index)
print len(samples)


# -- Survival
# px, py = 'COG3', 'COG2'
def survival(px, py):
    y = proteomics.ix[[px, py]].mean().dropna()

    # samples_up = set(y[y > np.percentile(y, 75)].index)
    # samples_down = set(y[y < np.percentile(y, 25)].index)

    samples_up = set(y[y > 2].index)
    samples_down = set(y[y < -2].index)

    print px, py, len(samples_down), len(samples_up)

    if len(samples_down) > 0 and len(samples_up) > 0:
        logrank = logrank_test(
            clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'],
            clinical.ix[samples_up, 'VITAL_STATUS'], clinical.ix[samples_down, 'VITAL_STATUS']
        )

        res = {
            'px': px, 'py': py,
            'ttest': logrank.test_statistic, 'logrank': logrank.p_value,
        }
        print res

        return res

res = [survival(px, py) for px, py in ppairs[['px', 'py']].values]
res = DataFrame([i for i in res if i])
res['fdr'] = multipletests(res['logrank'], method='fdr_bh')[1]
print res.sort('logrank')

#
plot_df = res[res['logrank'] < .05]

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs, pos = GridSpec(len(plot_df), 1, hspace=.5), 0
# px, py = 'SUPT16H', 'SSRP1'
for px, py in plot_df[['px', 'py']].values:
    y = proteomics.ix[[px, py]].mean().dropna()

    samples_up = set(y[y > 2].index)
    samples_down = set(y[y < -2].index)
    samples_bkg = set(y.index).difference(samples_up.union(samples_down))

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

    ax.set_title('Px (%s) ~ Py (%s)\nLogrank test: %.2e' % (px, py, logrank.p_value))
    ax.set_xlabel('Timeline (days)')
    ax.set_ylabel('Survival fraction')

    pos += 1

plt.gcf().set_size_inches(4, 3 * len(plot_df))
plt.savefig('%s/reports/survival_protein_pairs.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# --
# c = 4
def survival_complex(c):
    y = p_complex.ix[c].abs()

    samples_up = set(y[y > 2].index)
    samples_down = set(y[y < 2].index)

    print px, py, len(samples_down), len(samples_up)

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

res_c = [survival_complex(c) for c in p_complex.index]
res_c = DataFrame([i for i in res_c if i])
res_c['fdr'] = multipletests(res_c['logrank'], method='fdr_bh')[1]
print res_c.sort('logrank')

#
plot_df = res_c[(res_c['logrank'] < .05) & ([len(corum[i].intersection(ppairs_signif)) > 0 for i in res_c['complex']])]

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs, pos = GridSpec(len(plot_df), 1, hspace=.5), 0
# c = 5176
for c in plot_df['complex']:
    y = p_complex.ix[c].dropna().abs()

    samples_up = set(y[y > 2].index)
    samples_down = set(y[y < 2].index)

    logrank = logrank_test(
        clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'],
        clinical.ix[samples_up, 'VITAL_STATUS'], clinical.ix[samples_down, 'VITAL_STATUS']
    )

    # Plot
    ax = plt.subplot(gs[pos])

    kmf = KaplanMeierFitter()

    kmf.fit(clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_up, 'VITAL_STATUS'], label='Up')
    kmf.plot(ax=ax, ci_force_lines=False, color=palette_cnv_number[1])

    kmf.fit(clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_down, 'VITAL_STATUS'], label='Background')
    kmf.plot(ax=ax, ci_force_lines=False, color=palette_cnv_number[0])

    sns.despine(ax=ax)

    ax.set_title('%s\nLogrank test: %.2e' % (corum_n[c], logrank.p_value))
    ax.set_xlabel('Timeline (days)')
    ax.set_ylabel('Survival fraction')

    pos += 1

plt.gcf().set_size_inches(4, 3 * len(plot_df))
plt.savefig('%s/reports/survival_protein_complexes.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
