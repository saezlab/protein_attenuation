#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from protein_attenuation import wd, palette_cnv_number, cnv_names, default_color, palette_cnv
from sklearn.metrics import roc_curve, auc
from protein_attenuation.utils import randomise_matrix, gkn
from pandas import DataFrame, Series, read_csv, concat, merge
from sklearn.cross_validation import StratifiedShuffleSplit
from sklearn.linear_model import SGDClassifier


# -- Import data-sets
# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
cnv = cnv.applymap(lambda x: 0 if -1 <= x <= 1 else x)
cnv = cnv.loc[:, (cnv != 0).sum() != 0]
print cnv

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()
print clinical

# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print transcriptomics

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print proteomics


# -- Overlap
genes = set(cnv.index).intersection(proteomics.index).intersection(transcriptomics.index)
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
print len(genes), len(samples)


# -- Covariates
samplesheet = Series.from_csv('%s/data/samplesheet.csv' % wd)
samplesheet = {k[:15]: v for k, v in samplesheet.to_dict().items()}

design = Series([samplesheet[i[:15]] for i in samples], index=samples)
design = design.str.get_dummies()

design = concat([design, Series({i: clinical_gender[i] for i in design.index}).str.get_dummies()], axis=1)
design['age'] = [clinical_age[i] for i in design.index]
design_dict = design.T.to_dict()
print list(design)


# -- Classifications
def classify_cnv(c):
    # Build dependent and independent variables
    y = cnv.ix[genes, samples].unstack()

    y[y != c] = 0
    y[y == c] = 1

    x = DataFrame({
        'transcriptomics': transcriptomics.ix[genes, samples].unstack(),
        'proteomics': proteomics.ix[genes, samples].unstack()
    })
    x = concat([x, design.ix[[i[0] for i in x.index]].set_index(x.index)], axis=1).dropna()

    y = y.ix[x.index]

    # Build classification models
    lms = []
    for train, test in StratifiedShuffleSplit(y, test_size=.2, n_iter=15):
        # Fit model
        lm = SGDClassifier().fit(x.ix[train], y[train])

        # Predict test
        pred = concat([DataFrame(lm.decision_function(x.ix[test]), columns=['score'], index=x.ix[test].index), y.ix[test].rename('TP')], axis=1)

        # Evaluate fitted model
        curve_fpr, curve_tpr, _ = roc_curve(pred['TP'], pred['score'])
        curve_auc = auc(curve_fpr, curve_tpr)

        print c, curve_auc

        lms.append((lm, dict(zip(*(list(x), lm.coef_[0]))), curve_fpr, curve_tpr, curve_auc))

    return lms

l_models = {c: classify_cnv(c) for c in [-2, 2]}
l_models = {c: DataFrame(l_models[c], columns=['lm', 'features', 'fpr', 'tpr', 'auc'])  for c in [-2, 2]}
print '[INFO] Done'


# -- ROC curves
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})

# Plot all roc curves
for c in l_models:
    for fpr, tpr in l_models[c][['fpr', 'tpr']].values:
        plt.plot(fpr, tpr, c=sns.light_palette(palette_cnv_number[c], reverse=True).as_hex()[2], lw=.1)

# Plot median roc curve
for c in l_models:
    fpr, tpr = l_models[c][l_models[c]['auc'] == l_models[c]['auc'].median()][['fpr', 'tpr']].values[0]

    plt.plot(
        fpr, tpr, c=palette_cnv_number[c], lw=2.5, label='%s (AROC %.2f)' % (cnv_names[c], l_models[c]['auc'].median()),
        path_effects=[pe.Stroke(linewidth=5, foreground='white'), pe.Normal()]
    )


plt.plot([0, 1], [0, 1], 'k--', lw=.3)
sns.despine(trim=True)

plt.xlabel('False positive rate')
plt.ylabel('True positive rate')

plt.legend(loc='lower right')

plt.gcf().set_size_inches(5, 5)
plt.savefig('%s/reports/classification_aroc_cnv.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Coefficients
plot_df = DataFrame(
    [np.append(c, i) for c in l_models for i in DataFrame([i for i in l_models[c]['features']]).unstack().reset_index().values],
    columns=['type', 'feature', 'index', 'coef']
)

plot_df['type'] = plot_df['type'].replace(-2, 'depletion').replace(2, 'amplification')
plot_df['feature'] = [i.lower() for i in plot_df['feature']]

order = ['transcriptomics', 'proteomics', 'brca', 'coread', 'hgsc', 'female', 'male', 'age']

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.FacetGrid(plot_df, legend_out=True, sharex=False, sharey=True, size=2.5)
g.map(sns.boxplot, 'coef', 'feature', 'type', palette=palette_cnv, orient='h', linewidth=.3, sym='', order=order)
g.map(sns.stripplot, 'coef', 'feature', 'type', palette=palette_cnv, orient='h', jitter=True, edgecolor='white', linewidth=.3, split=True, size=3, order=order)
g.map(plt.axvline, x=0, ls='--', lw=.3, c='gray')
g.add_legend()
g.despine(trim=True)
g.set_axis_labels('Coefficients', 'Features')
g.set_titles(row_template='{row_name}')
g.fig.subplots_adjust(wspace=.1, hspace=.4)
plt.savefig('%s/reports/classification_features_cnv.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
