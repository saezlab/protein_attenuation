#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import itertools as it
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
# Clinical data
clinical = read_csv('./data/tcga_clinical.csv', index_col=0).dropna(subset=['time', 'status'])
clinical = clinical[(clinical['admin.disease_code'] == 'ov') & (clinical['time'] < (365 * 10))]
print 'clinical', clinical.shape

c_activities = read_csv('./tables/protein_complexes_proteomics_activities.csv')
# c_activities = c_activities[c_activities['fdr'] < .05]
c_activities = pivot_table(c_activities, index='complex', columns='sample', values='z', fill_value=np.nan)
print 'c_activities', c_activities.shape

uniprot = read_uniprot_genename()
corum = get_complexes_dict()
corum = {k: {uniprot[p][0] for p in corum[k] if p in uniprot} for k in corum}
corum_n = get_complexes_name()

# Overlap
samples = set(c_activities).intersection(clinical.index)
print 'samples', len(samples)

# Tumor suppressor
tms = read_csv('./tables/Cancer5000.oncodriveROLE.0.3-0.7.txt', sep='\t').append(read_csv('./tables/HCD.oncodriveROLE.0.3-0.7.txt', sep='\t'))

cgenes = {
    'suppressor': set(tms[tms['oncodriveROLE'] == 'Loss of function']['SYM']),
    'oncogene': set(tms[tms['oncodriveROLE'] == 'Activating']['SYM'])
}

cgenes_complexes = {
    'suppressor': {c for c in corum if len(corum[c].intersection(cgenes['suppressor'])) > 0},
    'oncogene': {c for c in corum if len(corum[c].intersection(cgenes['oncogene'])) > 0}
}
cgenes_complexes['suppressor'] = cgenes_complexes['suppressor'].difference(cgenes_complexes['oncogene'])
cgenes_complexes['oncogene'] = cgenes_complexes['oncogene'].difference(cgenes_complexes['suppressor'])
print 'cgenes_complexes', len(cgenes_complexes)


# -- Import regression results
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')

associations = {(px, py) for px, py in ppairs_cnv[['px', 'py']].values}
print len(associations)


# -- Logrank test
# px = 5450
surv = {}
for px in c_activities.index:
    df = c_activities.ix[px, samples].dropna()

    # if len(set(it.permutations(corum[px], 2)).intersection(associations)) > 0:
    if len(corum[px].intersection(cgenes['suppressor'])) > 0 or len(corum[px].intersection(cgenes['oncogene'])) > 0:
        samples_up = set(df[df > (df.mean() + df.std())].index)
        samples_dw = set(df[df < (df.mean() - df.std())].index)
        samples_bg = set(df.index).difference(samples_up.union(samples_dw))

        # samples_up = set(df[df > 0].index)
        # samples_dw = set(df[df < 0].index)
        # samples_bg = set(samples).difference(samples_up.union(samples_dw))

        if len(samples_up) > 10 and len(samples_dw) > 10:
            logrank = logrank_test(
                clinical.ix[samples_up, 'time'], clinical.ix[samples_dw, 'time'],
                clinical.ix[samples_up, 'status'], clinical.ix[samples_dw, 'status']
            )
            print logrank

            surv[px] = {'pval': logrank.p_value, 't': logrank.test_statistic, 'name': corum_n[px]}

surv = DataFrame(surv).T
surv['fdr'] = multipletests(surv['pval'], method='fdr_bh')[1]
print surv.sort(['fdr', 'pval'])

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
ax = plt.subplot(111)

kmf = KaplanMeierFitter()

kmf.fit(clinical.ix[samples_up, 'time'], event_observed=clinical.ix[samples_up, 'status'], label='Up (N=%d)' % len(samples_up))
kmf.plot(ci_force_lines=False, color=palette_cnv_number[1], ax=ax, show_censors=True, censor_styles={'ms': 5, 'marker': 'x'}, lw=1., ci_show=False)

kmf.fit(clinical.ix[samples_bg, 'time'], event_observed=clinical.ix[samples_bg, 'status'], label='Bkg (N=%d)' % len(samples_bg))
kmf.plot(ci_force_lines=False, color=palette_cnv_number[0], ax=ax, show_censors=True, censor_styles={'ms': 5, 'marker': 'x'}, lw=1., ci_show=False)

kmf.fit(clinical.ix[samples_dw, 'time'], event_observed=clinical.ix[samples_dw, 'status'], label='Down (N=%d)' % len(samples_dw))
kmf.plot(ci_force_lines=False, color=palette_cnv_number[-1], ax=ax, show_censors=True, censor_styles={'ms': 5, 'marker': 'x'}, lw=1., ci_show=False)

sns.despine(ax=ax)

ax.set_title('%s\nLogrank test: %.2e' % (corum_n[px], logrank.p_value))
ax.set_xlabel('Timeline (days)')
ax.set_ylabel('Survival fraction')

plt.gcf().set_size_inches(4, 2)
plt.savefig('./reports/survival_protein_regulators.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
