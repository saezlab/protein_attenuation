import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from sklearn.linear_model.ridge import Ridge
from scipy.stats import pearsonr, spearmanr
from pymist.utils.stringdb import get_stringdb
from pymist.utils.corumdb import get_complexes_pairs, get_complexes_dict
from pandas import DataFrame, Series, read_csv
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import p-p interactions
p_scores = read_csv('%s/tables/correlation_matrix.csv' % wd, index_col=0).replace(np.nan, 0)
print p_scores

# -- Import normalised pancan proteomics
pancan = read_csv('%s/tables/pancan_normalised.csv' % wd, index_col=0)
pancan = pancan[pancan.count(1) > (pancan.shape[1] * .25)]
print pancan


# -- Estimate complex activity
# c = 'SNARE complex (VAMP2, SNAP25, STX1a, STX3, CPLX1, CPLX3, CPLX4)'
# sample = 'TCGA-E2-A15A-01A-41-A21W-30'
def calc_activity(sample):
    y = pancan.ix[p_scores.index, sample].dropna()

    x = p_scores.ix[y.index].dropna()
    x = x.loc[:, (x != 0).sum() > 1]

    lm = Ridge(alpha=.01).fit(x, y)

    return dict(zip(*(x.columns, lm.coef_)))

p_activity = DataFrame({s: calc_activity(s) for s in pancan})
print p_activity.head()

# -- Export
p_activity.to_csv('%s/tables/protein_activity.csv' % wd)
print '[INFO] Exported'
