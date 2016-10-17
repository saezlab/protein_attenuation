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
from statsmodels.stats.weightstats import CompareMeans, DescrStatsW
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

# proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
# proteomics = proteomics[proteomics.count(1) > proteomics.shape[1] * .5]
proteomics_dict = {s: proteomics[s].dropna().to_dict() for s in proteomics}
print 'proteomics', proteomics.shape

# Overlap
samples, proteins = set(clinical.index).intersection(proteomics), set(proteomics.index)
print 'samples', 'proteins', len(samples), len(proteins)


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


# -- Protein complexes
uniprot = read_uniprot_genename()

corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if 1 < len(corum[k])}
corum_n = get_complexes_name()
print 'corum', len(corum)


# -- Estimate complex activity
def ztest_complex(s, c, df):
    x1 = [df[s][k] for k in df[s] if k in corum[c]]
    x2 = [df[s][k] for k in df[s] if k not in corum[c]]

    if len(x1) > 1:
        stat = CompareMeans(DescrStatsW(x1), DescrStatsW(x2))

        z_larger, p_larger = stat.ztest_ind(alternative='larger')
        z_smaller, p_smaller = stat.ztest_ind(alternative='smaller')

        z, p = z_larger, p_larger if p_larger < p_smaller else p_smaller

        res = {
            'sample': s, 'complex': c, 'name': corum_n[c],
            'z': z, 'pval': p, 'mean': np.mean(x1), 'targets': len(x1)
        }

        return res

c_activity = [ztest_complex(s, c, proteomics_dict) for s in samples for c in corum]
c_activity = DataFrame([i for i in c_activity if i])
c_activity['fdr'] = multipletests(c_activity['pval'],  method='fdr_bh')[1]
print c_activity.sort('fdr')

# Generate matrix from flaten list
c_activity_matrix = pivot_table(c_activity, index='complex', columns='sample', values='z', fill_value=np.nan)
print 'c_activity_matrix', c_activity_matrix.shape


# -- Logrank test
# px = 368
surv = {}
for px in c_activity_matrix.index:
    df = c_activity_matrix.ix[px, samples].dropna()

    if len(set(it.permutations(corum[px], 2)).intersection(associations)) > 0:
        samples_up = set(df[df > (df.mean() + df.std())].index)
        samples_dw = set(df[df < (df.mean() - df.std())].index)
        samples_bg = set(df.index).difference(samples_up.union(samples_dw))

        if len(samples_up) > len(samples) * .10 and len(samples_dw) > len(samples) * .10:
            logrank = logrank_test(
                clinical.ix[samples_up, 'time'], clinical.ix[samples_dw, 'time'],
                clinical.ix[samples_up, 'status'], clinical.ix[samples_dw, 'status']
            )
            print logrank

            surv[px] = {
                'pval': logrank.p_value, 't': logrank.test_statistic, 'name': corum_n[px], 'len': len(corum[px].intersection(proteins))
            }

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
