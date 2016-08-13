import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.preprocessing.imputation import Imputer
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pandas import DataFrame, Series, read_csv, pivot_table, concat

# -- Imports
corum_n = get_complexes_name()

# Proteomics complex activities
c_prot = read_csv('%s/tables/protein_complexes_proteomics_activities.tsv' % wd, sep='\t')
c_prot = pivot_table(c_prot, index='complex', columns='sample', values='z')
print 'c_prot', c_prot.shape

# Transcriptomics complex activities
c_trans = read_csv('%s/tables/protein_complexes_transcriptomics_activities.tsv' % wd, sep='\t')
c_trans = pivot_table(c_trans, index='complex', columns='sample', values='z')
print 'c_trans', c_trans.shape

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0).dropna(subset=['VITAL_STATUS', 'DAYS_TO_LAST_FOLLOWUP'])
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
clinical = clinical[clinical['DAYS_TO_LAST_FOLLOWUP'] > 1]
clinical = clinical.T
print 'clinical', clinical.shape


# -- Overlap
samples = set(c_prot).intersection(c_trans).intersection(clinical)
complexes = set(c_prot.index).intersection(c_trans.index)
c_prot, c_trans, clinical = c_prot.ix[complexes, samples].T, c_trans.ix[complexes, samples].T, clinical.ix[:, samples].T
print 'c_prot', 'c_trans', 'clinical', c_prot.shape, c_trans.shape, clinical.shape


# -- Impute missing values with mean
c_prot = DataFrame(Imputer().fit_transform(c_prot), index=c_prot.index, columns=c_prot.columns)

# -- Append clinical information
c_prot = concat([c_prot, clinical[['VITAL_STATUS', 'DAYS_TO_LAST_FOLLOWUP']]], axis=1)
c_trans = concat([c_trans, clinical[['VITAL_STATUS', 'DAYS_TO_LAST_FOLLOWUP']]], axis=1)


# -- Hazard proportional regression
# c = 285
def complex_survival(c):
    t_cf = CoxPHFitter(normalize=False)
    t_cf.fit(c_trans[[c, 'DAYS_TO_LAST_FOLLOWUP', 'VITAL_STATUS']], 'DAYS_TO_LAST_FOLLOWUP', event_col='VITAL_STATUS')
    t_cindex = concordance_index(c_trans['DAYS_TO_LAST_FOLLOWUP'], -t_cf.predict_partial_hazard(c_trans[[c]]).values.ravel(), c_trans['VITAL_STATUS'])
    # print t_cf.print_summary()

    p_cf = CoxPHFitter(normalize=False)
    p_cf.fit(c_prot[[c, 'DAYS_TO_LAST_FOLLOWUP', 'VITAL_STATUS']], 'DAYS_TO_LAST_FOLLOWUP', event_col='VITAL_STATUS')
    p_cindex = concordance_index(c_prot['DAYS_TO_LAST_FOLLOWUP'], -t_cf.predict_partial_hazard(c_prot[[c]]).values.ravel(), c_prot['VITAL_STATUS'])
    # print p_cf.print_summary()

    res = {'complex': c, 'p_cindex': p_cindex, 't_cindex': t_cindex}
    print res

    return res

c_survival = DataFrame([complex_survival(c) for c in complexes])
print c_survival.sort('p_cindex')


# -- Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})

plt.scatter(c_survival['t_cindex'], c_survival['p_cindex'])
plt.plot([0, 1], [0, 1], 'k--', lw=.3)
# sns.despine(trim=True)
plt.legend(loc='lower right')
plt.ylabel('Proteomics ~ Survival')
plt.xlabel('Transcriptomics ~ Survival')
plt.xlim(.4, .7)
plt.ylim(.4, .7)
plt.gcf().set_size_inches(4, 4)
plt.savefig('%s/reports/complex_survival_regression_scatter.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
