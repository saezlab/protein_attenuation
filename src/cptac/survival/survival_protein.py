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
from pandas import DataFrame, Series, read_csv, concat, merge


# -- Imports
# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
print cnv

# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print transcriptomics

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print proteomics

# Residuals
residuals = read_csv('%s/tables/protein_residuals.csv' % wd, index_col=0)
print residuals

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t', index_col=0).dropna(subset=['VITAL_STATUS', 'DAYS_TO_LAST_FOLLOWUP'])
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
clinical = clinical[clinical['DAYS_TO_LAST_FOLLOWUP'] > 1]
print clinical


# -- Overlap
genes = set(cnv.index).intersection(proteomics.index).intersection(transcriptomics.index).intersection(residuals.index)
samples = set(cnv).intersection(proteomics).intersection(transcriptomics).intersection(residuals).intersection(clinical.index)
print len(genes), len(samples)


# -- Survival
# c, p = 2, 'ERBB2'
def survival(p, c, thres=5):
    y = cnv.ix[p, samples].dropna()

    if len(y[y == c]) >= thres:
        logrank = logrank_test(
            clinical.ix[y[y == c].index, 'DAYS_TO_LAST_FOLLOWUP'], clinical.ix[y[y != c].index, 'DAYS_TO_LAST_FOLLOWUP'],
            clinical.ix[y[y == c].index, 'VITAL_STATUS'], clinical.ix[y[y != c].index, 'VITAL_STATUS']
        )

        print p, c, logrank.p_value

        return p, c, logrank.p_value
    return np.nan, np.nan, np.nan

res = DataFrame([survival(p, c) for p in genes for c in [-2, 2]], columns=['protein', 'type', 'pval']).dropna()
res['fdr'] = multipletests(res['pval'], method='fdr_bh')[1]
print res.sort('fdr')

# c = 2
# p = 'PCCB'
# y = cnv.ix[p, samples].dropna()
#
# ax = plt.subplot(111)
#
# kmf = KaplanMeierFitter()
#
# kmf.fit(clinical.ix[y[y == c].index, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[y[y == c].index, 'VITAL_STATUS'], label='%.0f' % c)
# kmf.plot(ax=ax, ci_force_lines=True)
#
# kmf.fit(clinical.ix[y[y != c].index, 'DAYS_TO_LAST_FOLLOWUP'], event_observed=clinical.ix[y[y != c].index, 'VITAL_STATUS'], label='%.0f' % .0)
# kmf.plot(ax=ax, ci_force_lines=True)


# -- Significant associations
a_transcriptomics = read_csv('%s/tables/regressions_overlap_transcriptomics_cnv.csv' % wd)
a_proteomics = read_csv('%s/tables/regressions_overlap_proteomics_cnv.csv' % wd)
a_residuals = read_csv('%s/tables/regressions_overlap_residuals_cnv.csv' % wd)

associations = {
    n: {(p, t) for p, t, e, f in df[['protein', 'type', 'effectsize', 'fdr']].values if f < 1e-2 and abs(e) > 1.}
    for n, df in [('Transcriptomics', a_transcriptomics), ('Proteomics', a_proteomics), ('Residuals', a_residuals)]
}
associations['Overlap'] = associations['Transcriptomics'].intersection(associations['Proteomics']).difference(associations['Residuals'])

for a in associations:
    res[a] = [int((p, t) in associations[a]) for p, t in res[['protein', 'type']].values]
print res.sort('fdr')


sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})

for a in associations:
    curve_fpr, curve_tpr, _ = roc_curve(res[a], 1 - res['fdr'])
    plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (a, auc(curve_fpr, curve_tpr)), c=palette[a])

plt.plot([0, 1], [0, 1], 'k--', lw=.3)
sns.despine(trim=True)
plt.legend(loc='lower right')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.gcf().set_size_inches(3, 3)
plt.savefig('%s/reports/survival_aroc.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
