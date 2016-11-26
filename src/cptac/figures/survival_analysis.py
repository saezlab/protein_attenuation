#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
from cptac import palette
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from pandas import read_csv


# -- Samples classification
cors = read_csv('./tables/samples_correlations.csv', index_col=0)

# -- Imports
# Clinical data
clinical = read_csv('./data/tcga_clinical.csv', index_col=0).dropna(subset=['time', 'status'])
clinical = clinical[(clinical['admin.disease_code'] == 'ov') & (clinical['time'] < (365 * 15))]
clinical = clinical[[i in cors.index for i in clinical.index]]
clinical['cluster'] = [cors.ix[i, 'cluster'] for i in clinical.index]
print 'clinical', clinical.shape


# -- Logrank test
logrank = logrank_test(
    clinical.ix[clinical['cluster'] == 1, 'time'], clinical.ix[clinical['cluster'] == 0, 'time'],
    clinical.ix[clinical['cluster'] == 1, 'status'], clinical.ix[clinical['cluster'] == 0, 'status']
)
print logrank

# -- Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
ax = plt.subplot(111)

kmf = KaplanMeierFitter()

kmf.fit(clinical.ix[clinical['cluster'] == 1, 'time'], event_observed=clinical.ix[clinical['cluster'] == 1, 'status'], label='Attenuated (N=%d)' % len(clinical.ix[clinical['cluster'] == 1, 'time']))
kmf.plot(ci_force_lines=False, color=palette['Transcriptomics'], ax=ax, show_censors=True, censor_styles={'ms': 5, 'marker': 'x'}, lw=1., ci_show=False)

kmf.fit(clinical.ix[clinical['cluster'] == 0, 'time'], event_observed=clinical.ix[clinical['cluster'] == 0, 'status'], label='Not attenuated (N=%d)' % len(clinical.ix[clinical['cluster'] == 0, 'time']))
kmf.plot(ci_force_lines=False, color=palette['Clinical'], ax=ax, show_censors=True, censor_styles={'ms': 5, 'marker': 'x'}, lw=1., ci_show=False)

sns.despine(ax=ax)

ax.set_title('Logrank test p-value: %.2e' % logrank.p_value)
ax.set_xlabel('Timeline (days)')
ax.set_ylabel('Survival fraction')

plt.gcf().set_size_inches(4, 2)
plt.savefig('./reports/survival_protein_regulators.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/survival_protein_regulators.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
