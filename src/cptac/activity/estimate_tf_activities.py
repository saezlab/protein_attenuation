import numpy as np
from cptac import wd
from sklearn.linear_model import Ridge
from pandas import DataFrame, read_csv


# -- TFs regulons
regulons = read_csv('%s/files/tfs_regulons.csv' % wd)
regulons = regulons.groupby('set')['gene'].agg(lambda x: set(x)).to_dict()
regulons = {s: {g: 1 for g in regulons[s]} for s in regulons}
regulons = DataFrame(regulons).replace(np.nan, 0).astype(np.int)
print regulons


# -- Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0).ix[regulons.index].dropna()
print transcriptomics


# -- Estimate activity: linear regression
# s = 'TCGA-E2-A15A-01'
tf_lms = {}
for s in transcriptomics.columns:
    y = transcriptomics[s].dropna()

    x = regulons.ix[y.index].dropna().astype(int)
    x = x.loc[:, x.sum() > 0]

    y = y.ix[x.index]

    lm = Ridge(alpha=.01).fit(x, y)
    print lm

    tf_lms[s] = dict(zip(*(x.columns, lm.coef_)))


# -- Export
tf_activities = DataFrame({tf: tf_lms[tf] for tf in tf_lms})
tf_activities.to_csv('%s/tables/tf_activities.csv' % wd)
print tf_activities
