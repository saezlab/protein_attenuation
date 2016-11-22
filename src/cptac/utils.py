#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import itertools as it
from cptac import wd
from scipy import stats
from pandas import read_csv, Series
from scipy.stats.distributions import hypergeom


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


def log_likelihood(y_true, y_pred):
    n = len(y_true)
    ssr = np.power(y_true - y_pred, 2).sum()
    var = ssr / n

    l = np.longfloat(1 / (np.sqrt(2 * np.pi * var))) ** n * np.exp(-(np.power(y_true - y_pred, 2) / (2 * var)).sum())
    ln_l = np.log(l)

    return ln_l


def f_statistic(y_true, y_pred, n, p):
    msm = np.power(y_pred - y_true.mean(), 2).sum() / p
    mse = np.power(y_true - y_pred, 2).sum() / (n - p - 1)

    f = msm / mse

    f_pval = stats.f.sf(f, p, n - p - 1)

    return f, f_pval


def r_squared(y_true, y_pred):
    sse = np.power(y_true - y_pred, 2).sum()
    sst = np.power(y_true - y_true.mean(), 2).sum()

    r = 1 - sse / sst

    return r


def ztest(targets, mu, var):
    z = (np.mean(targets) - mu) / (np.sqrt(var / len(targets)))
    p = 2 * stats.norm.sf(abs(z))
    return z, p, np.mean(targets), len(targets)


def read_gmt(file_path):
    with open(file_path) as f:
        signatures = {l.split('\t')[0]: set(l.strip().split('\t')[2:]) for l in f.readlines()}

    return signatures


def hypergeom_test(signature, background, sublist):
    """
    Performs hypergeometric test

    Arguements:
            signature: {string} - Signature IDs
            background: {string} - Background IDs
            sublist: {string} - Sub-set IDs

    # hypergeom.sf(x, M, n, N, loc=0)
    # M: total number of objects,
    # n: total number of type I objects
    # N: total number of type I objects drawn without replacement

    """
    pvalue = hypergeom.sf(
        len(sublist.intersection(signature)),
        len(background),
        len(background.intersection(signature)),
        len(sublist)
    )

    intersection = len(sublist.intersection(signature).intersection(signature))

    return pvalue, intersection


def jaccard(x, y):
    return float(len(x.intersection(y))) / len(x.union(y))


def import_signor():
    tab = read_csv('./tables/human_phosphorylations_26_09_16.csv', sep=';')
    tab = {(a, b) for a, b in tab[['ENTITYA', 'ENTITYB']].values}
    return tab


def import_kegg():
    tab = read_gmt('./tables/c2.cp.kegg.v5.1.symbols_only_metabolism.gmt')
    tab = {(a, b) for p in tab for a, b in it.combinations(tab[p], 2)}
    return tab