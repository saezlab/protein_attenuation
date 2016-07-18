import numpy as np
from cptac import wd
from scipy import stats
from pandas import read_csv, Series


def get_cancer_genes():
    types = ['Loss of function', 'Activating']

    table1 = read_csv('%s/files/Cancer5000.oncodriveROLE.0.3-0.7.txt' % wd, sep='\t')
    table2 = read_csv('%s/files/HCD.oncodriveROLE.0.3-0.7.txt' % wd, sep='\t')

    cancer_genes = {t: set(table1[table1['oncodriveROLE'] == t]['SYM']).union(table2[table2['oncodriveROLE'] == t]['SYM']) for t in types}

    return cancer_genes


def gkn(values):
    kernel = stats.gaussian_kde(values)
    return Series({k: np.log(kernel.integrate_box_1d(-1e4, v) / kernel.integrate_box_1d(v, 1e4)) for k, v in values.to_dict().items()})


def randomise_matrix(matrix):
    random_df = matrix.copy()
    movers = ~np.isnan(random_df.values)
    random_df.values[movers] = np.random.permutation(random_df.values[movers])
    return random_df
