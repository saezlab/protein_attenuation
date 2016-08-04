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
from pymist.utils.corumdb import get_complexes_name
from pandas import DataFrame, Series, read_csv, concat, merge, pivot_table


# -- Imports
# Protein complexes activities
p_complex_list = read_csv('%s/tables/protein_complexes_activities.tsv' % wd, sep='\t', index_col=0)
p_complex_z = pivot_table(p_complex_list, index='complex', columns='sample', values='z')
p_complex_n = get_complexes_name()
print p_complex_list.sort('fdr')

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0).dropna(subset=['VITAL_STATUS', 'DAYS_TO_LAST_FOLLOWUP'])
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
clinical = clinical[clinical['DAYS_TO_LAST_FOLLOWUP'] > 1]
print clinical


# -- Overlap
samples = set(p_complex_z).intersection(clinical.index)
print len(samples)


# -- Survival
# p_complex = 4
def survival(p_complex):
    samples_down = set(p_complex_list.ix[(p_complex_list['complex'] == p_complex) & (p_complex_list['fdr'] < .1) & (p_complex_list['z'] < 0), 'sample'])
    samples_up = set(p_complex_list.ix[(p_complex_list['complex'] == p_complex) & (p_complex_list['fdr'] < .1) & (p_complex_list['z'] > 0), 'sample'])
    print p_complex, len(samples_down), len(samples_up)

    if len(samples_down) > 0 and len(samples_up) > 0:
        logrank = logrank_test(
            clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'],
            clinical.ix[samples_up, 'VITAL_STATUS'], clinical.ix[samples_down, 'VITAL_STATUS']
        )

        res = {'complex': p_complex, 't_test': logrank.test_statistic, 'logrank': logrank.p_value}
        print res

        return res

res = [survival(p_complex) for p_complex in p_complex_z.index]
res = DataFrame([i for i in res if i])
res['fdr'] = multipletests(res['logrank'], method='fdr_bh')[1]
res['name'] = [p_complex_n[i] for i in res['complex']]
print res.sort('logrank')


# --
y = df.ix[res.ix[res['logrank'].argmin(), 'var'], samples].dropna()

samples_up = set(y[y > np.percentile(y, 75)].index)
samples_down = set(y[y < np.percentile(y, 25)].index)

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
ax = plt.subplot(111)

kmf = KaplanMeierFitter()

kmf.fit(clinical.ix[samples_up, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_up, 'VITAL_STATUS'], label='UP')
kmf.plot(ax=ax, ci_force_lines=False)

kmf.fit(clinical.ix[samples_down, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[samples_down, 'VITAL_STATUS'], label='DOWN')
kmf.plot(ax=ax, ci_force_lines=False)

sns.despine()

plt.gcf().set_size_inches(5, 3)
plt.savefig('%s/reports/survival_complexes.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
