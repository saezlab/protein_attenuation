import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, palette_cnv_number, cnv_names
from sklearn.metrics import roc_curve, auc
from cptac.utils import randomise_matrix, gkn
from pandas import DataFrame, Series, read_csv, concat, merge
from sklearn.cross_validation import StratifiedShuffleSplit
from sklearn.linear_model import SGDClassifier


# -- CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
cnv.columns = [i[:15] for i in cnv]
cnv = cnv.applymap(lambda x: 0 if -1 <= x <= 1 else x)
cnv = cnv.loc[:, (cnv != 0).sum() != 0]
cnv_dict = cnv.to_dict()
print cnv


# -- Protein residuals
residuals = read_csv('%s/tables/protein_residuals.csv' % wd, index_col=0)
print residuals


# -- Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()


# -- Covariates
samplesheet = read_csv('%s/data/samplesheet.csv' % wd, index_col=0)
samplesheet = {k[:15]: v for k, v in samplesheet['type'].to_dict().items()}

design = Series([samplesheet[i[:15]] for i in residuals], index=residuals.columns)
design = design.str.get_dummies()

design = concat([design, Series({i: clinical_gender[i] for i in design.index}).str.get_dummies()], axis=1)
design['age'] = [clinical_age[i] for i in design.index]
# design['tmt'] = np.bitwise_or(design['BRCA'], design['HGSC'])
# design['shotgun'] = design['COREAD']
design_dict = design.T.to_dict()
print list(design)


# --
genes = set(cnv.index).intersection(residuals.index)
samples = set(cnv).intersection(residuals)

y = cnv.ix[genes, samples].unstack()

x = residuals.ix[genes, samples].unstack().dropna()
x.name = 'residuals'

x = concat([x, DataFrame({i: design_dict[i[0]] for i in x.index}).T], axis=1)

y = y.ix[x.index]

# [(train, test) for train, test in StratifiedShuffleSplit(y, test_size=.1)]
lms = [(train, test, SGDClassifier().fit(x.ix[train], y[train])) for train, test in StratifiedShuffleSplit(y, test_size=.1)]
print '[INFO] Done'


# -- Generate roc curves
classes_roc = []
for c in cnv_names:
    for train, test, lm in lms:
        predicted = concat([DataFrame(lm.decision_function(x.ix[test]), columns=lm.classes_, index=x.ix[test].index), y.ix[test].rename('TP')], axis=1)

        curve_fpr, curve_tpr, _ = roc_curve((predicted['TP'] == c).astype(int), predicted[c])
        curve_auc = auc(curve_fpr, curve_tpr)

        print c, curve_auc

        classes_roc.append((c, curve_fpr, curve_tpr, curve_auc))

classes_roc = DataFrame(classes_roc, columns=['class', 'fpr', 'tpr', 'auc'])
print classes_roc


# -- ROC curves
cnv_arocs = {c: classes_roc.loc[classes_roc['class'] == c, 'auc'].mean() for c in cnv_names}

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})

for c, fpr, tpr, auc in classes_roc.values:
    plt.plot(fpr, tpr, c=palette_cnv_number[c], label='%s (AUC %0.2f)' % (cnv_names[c], cnv_arocs[c]))

plt.plot([0, 1], [0, 1], 'k--', lw=.3)
sns.despine(trim=True)

plt.xlabel('False positive rate')
plt.ylabel('True positive rate')

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='lower right')

plt.gcf().set_size_inches(5, 5)
plt.savefig('%s/reports/pancan_residuals_cnv_aroc.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
