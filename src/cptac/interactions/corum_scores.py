import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pymist.utils.corumdb import get_complexes_dict
from pandas import DataFrame, Series, read_csv
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Uniprot
uniprot = read_uniprot_genename()


# -- Import CORUM
corum = get_complexes_dict()
corum = {k: {uniprot[i][0] for i in v if i in uniprot} for k, v in corum.items()}


# -- Estimate complex
def calc_activity(df, line):
    y = df.ix[k_targets.index, line].dropna()

    x = k_targets.ix[y.index].dropna()
    x = x.loc[:, x.sum() != 0]

    lm = Ridge(alpha=.01).fit(x, y)

    return dict(zip(*(x.columns, lm.coef_)))

k_activity = DataFrame({l: calc_activity(pproteomics, l) for l in pproteomics})
print k_activity.head()