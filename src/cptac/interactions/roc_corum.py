import itertools as it
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, palette
from gdsc import wd as gdsc_wd
from mtkirc.utils import gkn
from pandas.stats.misc import zscore
from sklearn.linear_model.base import LinearRegression
from pandas import read_csv, DataFrame, concat, Series
from scipy.stats.stats import spearmanr
from statsmodels.distributions import ECDF
from sklearn.metrics.ranking import roc_curve, auc
from pymist.utils.corumdb import get_complexes_pairs
from statsmodels.stats.multitest import multipletests
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Uniprot
uniprot = read_uniprot_genename()


# -- Import CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
proteins = {i for p1, p2 in corum for i in [p1, p2]}


# -- Imports
# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()

# Proteomics
brca = read_csv('%s/tables/brca_proteomics_processed.csv' % wd, index_col=0)
brca = brca[brca.count(1) > (brca.shape[1] * .5)]

hgsc = read_csv('%s/tables/hgsc_proteomics_processed.csv' % wd, index_col=0)
hgsc = hgsc[hgsc.count(1) > (hgsc.shape[1] * .5)]

coread = read_csv('%s/tables/coread_proteomics_processed.csv' % wd, index_col=0)
coread = coread[coread.count(1) > (coread.shape[1] * .5)]

ov_prot = set(brca.index).intersection(hgsc.index).intersection(coread.index).intersection(proteins)
brca, hgsc, coread = brca.ix[ov_prot], hgsc.ix[ov_prot], coread.ix[ov_prot]

pancan = read_csv('%s/tables/pancan_preprocessed_normalised.csv' % wd, index_col=0).ix[ov_prot]


# -- Protein-protein correlation
cor_res = {}
for name, d_df in [('BRCA', brca), ('HGSC', hgsc), ('COREAD', coread), ('Pancancer', pancan)]:
    df = d_df.T.corr(method='pearson')
    df.values[np.tril_indices(df.shape[0], 0)] = np.nan
    df.index.name = None
    df = df.unstack().reset_index().dropna()
    df.columns = ['p1', 'p2', 'cor']

    df['tp'] = [1 if ((p1, p2) in corum) else 0 for p1, p2 in df[['p1', 'p2']].values]
    df['score'] = df['cor'].abs()

    print df

    curve_fpr, curve_tpr, _ = roc_curve(df['tp'], df['score'])

    cor_res[name] = (curve_fpr, curve_tpr)

    print '[INFO] Correlation: %d/%d' % (df['tp'].sum(), len(df))


# -- Plot
sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=0.75)
for name in ['BRCA', 'COREAD', 'HGSC', 'Pancancer']:
    curve_fpr, curve_tpr = cor_res[name]
    plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (name, auc(curve_fpr, curve_tpr)), c=palette[name])

plt.plot([0, 1], [0, 1], 'k--', lw=.3)
sns.despine(trim=True)
plt.legend(loc='lower right')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('CORUM - ROC curve')
plt.gcf().set_size_inches(3, 3)
plt.savefig('%s/reports/corum_auc.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
