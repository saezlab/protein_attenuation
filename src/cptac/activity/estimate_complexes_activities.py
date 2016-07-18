import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from cptac.utils import gkn
from sklearn.linear_model import Ridge, LinearRegression
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Complexes proteins
uniprot = read_uniprot_genename()
corum_dict = {k: {uniprot[p][0]: 1 for p in v if p in uniprot} for k, v in get_complexes_dict().items()}
corum_dict = {k: corum_dict[k] for k in corum_dict if len(corum_dict[k]) > 1}
corum = DataFrame(corum_dict).replace(np.nan, 0).astype(np.int)
corum_names = get_complexes_name()
print corum


# -- Proteomics
proteomics = read_csv('%s/data/pancan_preprocessed.csv' % wd, index_col=0)
proteomics.columns = [i[:15] for i in proteomics]

remove_samples = {i for i in set(proteomics) if proteomics.loc[:, [i]].shape[1] == 2 and proteomics.loc[:, [i]].corr().ix[0, 1] < .4}
proteomics = proteomics.drop(remove_samples, axis=1)

proteomics = DataFrame({i: proteomics.loc[:, [i]].mean(1) for i in set(proteomics)})
proteomics = proteomics[proteomics.count(1) > (proteomics.shape[1] * .5)]
print proteomics


# -- Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()


# -- Covariates
samplesheet = read_csv('%s/data/samplesheet.csv' % wd, index_col=0)
samplesheet = {k[:15]: v for k, v in samplesheet['type'].to_dict().items()}

design = Series([samplesheet[i[:15]] for i in proteomics], index=proteomics.columns)
design = design.str.get_dummies()

design = concat([design, Series({i: clinical_gender[i] for i in design.index}).str.get_dummies()], axis=1)
design['age'] = [clinical_age[i] for i in design.index]
design['tmt'] = np.bitwise_or(design['BRCA'], design['HGSC'])
design['shotgun'] = design['COREAD']
design_dict = design.T.to_dict()
print list(design)


# -- Normalise pancan data-set
def rm_batch(x, y, covariates=['BRCA', 'COREAD', 'HGSC', 'FEMALE', 'MALE', 'age', 'tmt', 'shotgun']):
    ys = y.dropna()
    xs = x.ix[ys.index, covariates]

    lm = LinearRegression().fit(xs, ys)

    return ys - xs.dot(lm.coef_) - lm.intercept_

proteomics = DataFrame({p: rm_batch(design, proteomics.ix[p, design.index]) for p in proteomics.index}).T
print '[INFO] Covariates regressed-out'


# -- Centre gene
proteomics = DataFrame({i: gkn(proteomics.ix[i].dropna()).to_dict() for i in proteomics.index}).T


# -- Estimate activity: linear regression
# s = 'TCGA-E2-A15A-01'
complexes = {}
for s in proteomics.columns:
    y = proteomics[s].dropna()

    x = corum.ix[y.index].dropna().astype(int)
    x = x.loc[:, x.sum() > 1]

    y = y.ix[x.index]

    lm = Ridge(alpha=.01).fit(x, y)
    print lm

    complexes[s] = dict(zip(*(x.columns, lm.coef_)))


# -- Export
complex_activities = DataFrame({tf: complexes[tf] for tf in complexes})
complex_activities.to_csv('%s/tables/protein_complexes_activities.csv' % wd)
print complex_activities
