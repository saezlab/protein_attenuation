import itertools as it
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from gdsc import wd as gdsc_wd
from mtkirc.utils import gkn
from pandas.stats.misc import zscore
from pandas import read_csv, DataFrame, concat
from scipy.stats.stats import spearmanr
from sklearn.metrics.ranking import roc_curve, auc
from pymist.utils.corumdb import get_complexes_pairs
from statsmodels.stats.multitest import multipletests
from pymist.utils.map_peptide_sequence import read_uniprot_genename

# -- Uniprot
uniprot = read_uniprot_genename()

# -- Import CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}

# -- Imports
brca = read_csv('%s/tables/brca_proteomics_processed.csv' % wd, index_col=0)
coread = read_csv('%s/tables/coread_proteomics_processed.csv' % wd, index_col=0)
hgsc = read_csv('%s/tables/hgsc_proteomics_processed.csv' % wd, index_col=0)

coread_c = read_csv('%s/data/proteomics_coread_processed.csv' % gdsc_wd, index_col=0)
coread_c = coread_c[[i in uniprot for i in coread_c.index]]
coread_c['id'] = [uniprot[i][0] for i in coread_c.index]
coread_c = coread_c.groupby('id').mean()

bcra_c = read_csv('%s/data/proteomics_brca_processed.csv' % gdsc_wd, index_col=0)
bcra_c = bcra_c[[i in uniprot for i in bcra_c.index]]
bcra_c['id'] = [uniprot[i][0] for i in bcra_c.index]
bcra_c = bcra_c.groupby('id').mean()

ov_prot = set(brca.index).intersection(coread.index).intersection(hgsc.index).intersection(coread_c.index).intersection(bcra_c.index)
brca, coread, hgsc, coread_c, bcra_c = brca.ix[ov_prot], coread.ix[ov_prot], hgsc.ix[ov_prot], coread_c.ix[ov_prot], bcra_c.ix[ov_prot]

all = concat([brca, coread, hgsc, coread_c, bcra_c], axis=1)
all_n = zscore(all)

# -- Protein-protein correlation
cor_res = {}
for name, d_df in [('BRCA', brca), ('COREAD', coread), ('HGSC', hgsc), ('All', all), ('All (normalised)', all_n), ('COREAD (lines)', coread_c), ('BRCA (lines)', bcra_c)]:
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
pal = dict(zip(*(['BRCA', 'COREAD', 'HGSC', 'All', 'All (normalised)', 'COREAD (lines)', 'BRCA (lines)'], sns.color_palette('Set2', 7).as_hex())))

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=0.75)
for name in ['All', 'All (normalised)', 'BRCA', 'COREAD', 'HGSC', 'COREAD (lines)', 'BRCA (lines)']:
    curve_fpr, curve_tpr = cor_res[name]
    plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (name, auc(curve_fpr, curve_tpr)), c=pal[name])

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
