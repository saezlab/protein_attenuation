import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from cptac import wd, palette
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
p_complex = read_csv('%s/tables/protein_complexes_activities_mean.tsv' % wd, sep='\t', index_col=0)
print p_complex

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0).dropna(subset=['VITAL_STATUS', 'DAYS_TO_LAST_FOLLOWUP'])
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
clinical = clinical[clinical['DAYS_TO_LAST_FOLLOWUP'] > 1]
print clinical

# P-pairs
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot} for k, v in get_complexes_dict().items()}
corum_n = get_complexes_name()

ppairs = read_csv('%s/tables/ppairs_cnv_regulation.csv' % wd)
ppairs_corum = {i: [(c, corum_n[c]) for c in corum if i in corum[c]] for i in ppairs['px']}
print ppairs


# -- Overlap
samples = set(p_complex).intersection(clinical.index)
print len(samples)


# -- Survival
# p = 4
def survival(p):
    samples_up = np.percentile(p_complex.ix[p].dropna(), 90)
    samples_down = np.percentile(p_complex.ix[p].dropna(), 10)

    samples_up = set(p_complex.ix[p][p_complex.ix[p] > samples_up].index)
    samples_down = set(p_complex.ix[p][p_complex.ix[p] < samples_down].index)
    samples_bkg = set(p_complex.ix[p].dropna().index).difference(samples_up.union(samples_down))
    print p, len(samples_down), len(samples_up), len(samples_bkg)

    if len(samples_down) > 0 and len(samples_up) > 0:
        logrank_up = logrank_test(
            clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[samples_bkg.union(samples_down), 'DAYS_TO_LAST_FOLLOWUP'],
            clinical.ix[samples_up, 'VITAL_STATUS'], clinical.ix[samples_bkg.union(samples_down), 'VITAL_STATUS']
        )

        logrank_down = logrank_test(
            clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[samples_bkg.union(samples_up), 'DAYS_TO_LAST_FOLLOWUP'],
            clinical.ix[samples_down, 'VITAL_STATUS'], clinical.ix[samples_bkg.union(samples_up), 'VITAL_STATUS']
        )

        res = {
            'complex': p,
            'up_ttest': logrank_up.test_statistic, 'up_logrank': logrank_up.p_value,
            'down_ttest': logrank_down.test_statistic, 'down_logrank': logrank_down.p_value
        }
        print res

        return res

res = [survival(p) for p in p_complex.index]
res = DataFrame([i for i in res if i])
res['up_fdr'] = multipletests(res['up_logrank'], method='fdr_bh')[1]
res['down_fdr'] = multipletests(res['down_logrank'], method='fdr_bh')[1]
res['name'] = [corum_n[i] for i in res['complex']]
print res.sort('up_fdr')


# --
p = 2383

samples_up = np.percentile(p_complex.ix[p].dropna(), 90)
samples_down = np.percentile(p_complex.ix[p].dropna(), 10)

samples_up = set(p_complex.ix[p][p_complex.ix[p] > samples_up].index)
samples_down = set(p_complex.ix[p][p_complex.ix[p] < samples_down].index)
samples_bkg = set(p_complex.ix[p].dropna().index).difference(samples_up.union(samples_down))
print p_complex, len(samples_down), len(samples_up), len(samples_bkg)

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
ax = plt.subplot(111)

kmf = KaplanMeierFitter()

kmf.fit(clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_up, 'VITAL_STATUS'], label='UP')
kmf.plot(ax=ax, ci_force_lines=False)

kmf.fit(clinical.ix[samples_bkg, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_bkg, 'VITAL_STATUS'], label='BKG')
kmf.plot(ax=ax, ci_force_lines=False)

kmf.fit(clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_down, 'VITAL_STATUS'], label='DOWN')
kmf.plot(ax=ax, ci_force_lines=False)

sns.despine()

plt.gcf().set_size_inches(5, 3)
plt.savefig('%s/reports/survival_protein_complexes.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
