import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from cptac.utils import gkn
from sklearn.linear_model import Ridge, LinearRegression
from pandas import DataFrame, Series, read_csv, concat


# -- TFs regulons
regulons = read_csv('%s/files/tfs_regulons.csv' % wd)
regulons = regulons.groupby('set')['gene'].agg(lambda x: set(x)).to_dict()
regulons = {s: {g: 1 for g in regulons[s]} for s in regulons}
regulons = DataFrame(regulons).replace(np.nan, 0).astype(np.int)
print regulons


# -- Gexp
transcriptomics = read_csv('%s/data/tcga_rnaseq.tsv' % wd, sep='\t', index_col=0).ix[regulons.index].dropna()
transcriptomics.columns = [i[:15] for i in transcriptomics]
transcriptomics = DataFrame({i: transcriptomics.loc[:, [i]].mean(1) for i in set(transcriptomics)})
print transcriptomics


# -- Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()


# -- Covariates
samplesheet = read_csv('%s/data/samplesheet.csv' % wd, index_col=0)
samplesheet = {k[:15]: v for k, v in samplesheet['type'].to_dict().items()}

design = Series([samplesheet[i[:15]] for i in transcriptomics], index=transcriptomics.columns)
design = design.str.get_dummies()

design = concat([design, Series({i: clinical_gender[i] for i in design.index}).str.get_dummies()], axis=1)
design['age'] = [clinical_age[i] for i in design.index]
design_dict = design.T.to_dict()
print list(design)


# -- Normalise pancan data-set
def rm_batch(x, y, covariates=['BRCA', 'COREAD', 'HGSC', 'FEMALE', 'MALE', 'age']):
    ys = y.dropna()
    xs = x.ix[ys.index, covariates]

    lm = LinearRegression().fit(xs, ys)

    return ys - xs.dot(lm.coef_) - lm.intercept_

transcriptomics = DataFrame({p: rm_batch(design, transcriptomics.ix[p, design.index]) for p in transcriptomics.index}).T
print '[INFO] Covariates regressed-out'


# -- Centre gene
transcriptomics = DataFrame({i: gkn(transcriptomics.ix[i].dropna()).to_dict() for i in transcriptomics.index}).T


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
