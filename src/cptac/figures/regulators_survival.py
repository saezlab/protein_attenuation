#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import palette_cnv_number
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index
from lifelines.statistics import logrank_test
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LinearRegression
from statsmodels.duration.hazard_regression import PHReg
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pandas import DataFrame, Series, read_csv, pivot_table, concat


# -- Imports
# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)
print 'cnv', cnv.shape

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape

# Clinical data
clinical = read_csv('./data/tcga_clinical.csv', index_col=0).dropna(subset=['time', 'status'])
clinical = clinical[(clinical['time'] > 1) & (clinical['time'] < 3650)]
clinical = clinical[clinical['admin.disease_code'] == 'ov']
print 'clinical', clinical.shape

# Overlap
genes = set(transcriptomics.index).intersection(proteomics.index)
samples = set(proteomics).intersection(cnv).intersection(clinical.index)
print 'samples', len(samples)


# -- Import regression results
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')


# -- Regressions residuals
# px = 'EIF3A'
res = {}
for px in set(set(ppairs_cnv['px'])):
    pys = list(set(ppairs_cnv.loc[ppairs_cnv['px'] == px, 'py']))

    df = concat([transcriptomics.ix[[px], samples].T, proteomics.ix[pys, samples].T], axis=1).dropna()

    if df.shape[0] > 0:
        lm = sm.OLS(df[px], sm.add_constant(df[pys])).fit()
        print lm.summary()

        res[px] = lm.resid.to_dict()

res = DataFrame(res).T
print res


# -- Survival analysis
# px = 'EEF1E1'
surv = {}
for px in res.index:
    df = res.ix[px, samples].dropna()

    samples_up = set(df[df > df.quantile(.75)].index)
    samples_dw = set(df[df < df.quantile(.25)].index)
    samples_bg = set(df.index).difference(samples_up.union(samples_dw))

    logrank = logrank_test(
        clinical.ix[samples_up, 'time'], clinical.ix[samples_dw, 'time'],
        clinical.ix[samples_up, 'status'], clinical.ix[samples_dw, 'status']
    )
    print logrank

    surv[px] = {'pval': logrank.p_value, 't': logrank.test_statistic}
surv = DataFrame(surv).T
surv['fdr'] = multipletests(surv['pval'], method='fdr_bh')[1]
print surv.sort(['fdr', 'pval'])


# Plot
sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'}, font_scale=0.75)
ax = plt.subplot(111)

kmf = KaplanMeierFitter()

kmf.fit(clinical.ix[samples_up, 'time'], event_observed=clinical.ix[samples_up, 'status'], label='Up')
kmf.plot(ci_force_lines=False, color=palette_cnv_number[1], ax=ax, show_censors=True, censor_styles={'ms': 5})

kmf.fit(clinical.ix[samples_bg, 'time'], event_observed=clinical.ix[samples_bg, 'status'], label='Bkg')
kmf.plot(ci_force_lines=False, color=palette_cnv_number[0], ax=ax, show_censors=True, censor_styles={'ms': 5})

kmf.fit(clinical.ix[samples_dw, 'time'], event_observed=clinical.ix[samples_dw, 'status'], label='Down')
kmf.plot(ci_force_lines=False, color=palette_cnv_number[-1], ax=ax, show_censors=True, censor_styles={'ms': 5})

sns.despine(ax=ax)

ax.set_title('Logrank test: %.2e' % logrank.p_value)
ax.set_xlabel('Timeline (days)')
ax.set_ylabel('Survival fraction')

plt.gcf().set_size_inches(5, 3)
plt.savefig('./reports/survival_protein_regulators.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
