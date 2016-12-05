#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from protein_attenuation import wd, palette_cnv_number, palette
from statsmodels.duration.hazard_regression import PHReg
from pandas import DataFrame, Series, read_csv, concat


# -- Import
# Proteomics samples
samples = set(Series.from_csv('%s/data/samplesheet.csv' % wd).index)
print len(samples)

# Clinical data
clinical = read_csv('%s/data/tcga_clinical.csv' % wd, index_col=0).dropna(subset=['time', 'status'])
clinical = clinical[clinical['time'] > 1]
clinical.index = ['%s-01' % i for i in clinical.index]
clinical = clinical.ix[samples].dropna(subset=['time', 'status'])
clinical['admin.disease_code'] = ['coread' if i in ['coad', 'read'] else i for i in clinical['admin.disease_code']]
clinical = concat([clinical, clinical['admin.disease_code'].str.get_dummies()], axis=1)
print 'clinical', clinical.shape


# --
plot_df = clinical[['brca', 'coread', 'ov', 'time', 'status']]

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'}, font_scale=0.75)
ax = plt.subplot(111)

kmf = KaplanMeierFitter()

kmf.fit(plot_df.ix[plot_df['brca'] == 1, 'time'], event_observed=plot_df.ix[plot_df['brca'] == 1, 'status'], label='BRCA')
kmf.plot(ci_force_lines=False, color=palette['BRCA'], ax=ax, show_censors=True, censor_styles={'ms': 5})

kmf.fit(plot_df.ix[plot_df['coread'] == 1, 'time'], event_observed=plot_df.ix[plot_df['coread'] == 1, 'status'], label='COREAD')
kmf.plot(ci_force_lines=False, color=palette['COREAD'], ax=ax, show_censors=True, censor_styles={'ms': 5})

kmf.fit(plot_df.ix[plot_df['ov'] == 1, 'time'], event_observed=plot_df.ix[plot_df['ov'] == 1, 'status'], label='HGSC')
kmf.plot(ci_force_lines=False, color=palette['HGSC'], ax=ax, show_censors=True, censor_styles={'ms': 5})

sns.despine(ax=ax)

ax.set_xlabel('Timeline (days)')
ax.set_ylabel('Survival fraction')

plt.gcf().set_size_inches(4, 3)
plt.savefig('%s/reports/survival_tissue_origin.pdf' % wd, bbox_inches='tight')
plt.savefig('%s/reports/survival_tissue_origin.png' % wd, bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'
