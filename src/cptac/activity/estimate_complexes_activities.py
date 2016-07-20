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
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print proteomics

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
