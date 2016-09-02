#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
from cptac import palette_cnv_number, palette
from cptac.utils import gkn, randomise_matrix
from pandas import DataFrame, Series, read_csv, concat
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index, k_fold_cross_validation
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests
from sklearn.cross_validation import ShuffleSplit, KFold
from sklearn.linear_model.base import LinearRegression
from statsmodels.duration.hazard_regression import PHReg
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict


# -- Clinical
clinical = read_csv('./data/hgsc_clinical.csv', index_col=0)
clinical.index = ['%s-01' % i for i in clinical.index]

clinical['time'] = [i for i in clinical['daystodeath_or_LFU']]
clinical['status'] = [0 if i == 'LIVING' else 1 for i in clinical['vital_status']]

clinical = clinical[clinical['time'] < (10 * 365)]
print 'clinical', clinical.shape


# -- Proteomics
proteomics = read_csv('./data/hgsc_proteomics_processed.csv', index_col=0)
proteomics.columns = [i[:15] for i in proteomics]

# Average replicates
remove_samples = {i for i in set(proteomics) if proteomics.loc[:, [i]].shape[1] == 2 and proteomics.loc[:, [i]].corr().ix[0, 1] < .4}
proteomics = proteomics.drop(remove_samples, axis=1)
proteomics = DataFrame({i: proteomics.loc[:, [i]].mean(1) for i in set(proteomics)})

# Drop missing values
proteomics = proteomics.dropna()

# Overlap
proteins, samples = list(set(proteomics.index)), list(set(proteomics).intersection(clinical.index))

# Normalise
proteomics = DataFrame({i: gkn(proteomics.ix[i].dropna()).to_dict() for i in proteins}).T


# -- Protein complexes
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if len(corum[k]) > 1}
corum_n = get_complexes_name()
print 'corum', len(corum)


# -- Signif Px -> Py associations
ppairs = read_csv('./tables/ppairs_cnv_regulation.csv')
p_regulators = {p for p in ppairs['px']}
p_complexes = list({(px, c) for px, py in ppairs[['px', 'py']].values for c in corum if px in corum[c] and py in corum[c]})
print 'p_complexes', len(p_complexes)


# -- Regression: c = 351
strata = concat([
    clinical['age_at_diagnosis'],
    clinical['tumor_stage'].str.get_dummies(),
    clinical['tumor_grade'].str.get_dummies(),
    clinical['PlatinumStatus'].str.get_dummies(),
    clinical['time'],
    clinical['status']
], axis=1)

m1 = CoxPHFitter(normalize=False, penalizer=.1).fit(strata, 'time', event_col='status', include_likelihood=True)
print m1.print_summary()
print m1._log_likelihood

c_surv = {}
for c in corum:
    x = concat([strata, proteomics.ix[corum[c], samples].T], axis=1).dropna()
    m2 = CoxPHFitter(normalize=False, penalizer=.1).fit(x, 'time', event_col='status', include_likelihood=True)

    print m2.print_summary()
    print m2._log_likelihood

    lr = 2 * (m2._log_likelihood - m1._log_likelihood)

    df = x.drop(['time', 'status'], axis=1).shape[1] - strata.drop(['time', 'status'], axis=1).shape[1]
    pval = stats.chi2.sf(lr, df)

    res = {
        'c': c, 'name': corum_n[c],
        'df': df, 'll_m1': m1._log_likelihood, 'll_m2': m2._log_likelihood,
        'lr': lr, 'pval': pval,
        'regulators': '_'.join(corum[c].intersection(p_regulators))
    }
    print pval, df

    c_surv[c] = res

c_surv = DataFrame(c_surv).T
c_surv['adj_pval'] = multipletests(c_surv['pval'], method='bonferroni')[1]
print c_surv.sort('pval')

print c_surv[[i != '' for i in c_surv['regulators']]].sort('pval')
