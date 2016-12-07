#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import itertools as it
from scipy import stats
from pandas import read_csv, Series
from scipy.stats.distributions import hypergeom


uniprot_fasta = './files/uniprot_sprot_2016_01_11.fasta'


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


# -- Uniprot related functions
def read_fasta(fasta_file=None, os='Homo sapiens'):
    fasta_file = uniprot_fasta if fasta_file is None else fasta_file

    sequences = {}

    with open(fasta_file) as f:
        lines = f.readlines()

        for i in range(len(lines)):
            if lines[i].startswith('>sp') and (('OS='+os) in lines[i]):
                uniprot = lines[i].split('|')[1].strip()
                sequence = ''

                i += 1
                while (i < len(lines)) and (not lines[i].startswith('>sp')):
                    sequence += lines[i].strip()
                    i += 1

                sequences[uniprot] = sequence

    return sequences


def match_sequence(sequences, sequence):
    return [k for k, v in sequences.items() if sequence in v]


def read_uniprot_genename(fasta_file=None, os='Homo sapiens'):
    fasta_file = uniprot_fasta if fasta_file is None else fasta_file

    uniprot2genename = {}

    with open(fasta_file) as f:
        lines = f.readlines()

        for i in range(len(lines)):
            if lines[i].startswith('>sp') and (('OS='+os) in lines[i]) and ('GN=' in lines[i]):
                uniprot = lines[i].split('|')[1].strip()
                genename = lines[i].split(' GN=')[1].split(' ')[0]
                accname = lines[i].split('|')[2].strip().split(' ')[0].strip()
                uniprot2genename[uniprot] = (genename, accname)

    return uniprot2genename


def read_uniprot_accname(fasta_file=None, os='Homo sapiens'):
    fasta_file = uniprot_fasta if fasta_file is None else fasta_file

    uniprot2accname = {}

    with open(fasta_file) as f:
        lines = f.readlines()

        for i in range(len(lines)):
            if lines[i].startswith('>sp') and (('OS='+os) in lines[i]):
                uniprot = lines[i].split('|')[1].strip()
                accname = lines[i].split('|')[2].strip().split(' ')[0].strip()
                uniprot2accname[uniprot] = accname

    return uniprot2accname


# -- CORUM complexes related functions
def get_complexes_pairs(organism='Human', corum_file='./files/allComplexesCore.csv'):
    complexes = read_csv(corum_file, sep=';').dropna(subset=['Complex id'])

    complexes = complexes[complexes['organism'] == organism]

    complexes = {c: {x.replace('(', '').replace(')', '') for p in complexes.loc[complexes['Complex id'] == c, 'subunits (UniProt IDs)'] for x in p.split(',')} for c in complexes['Complex id']}

    complexes_pairs = {(p1, p2) for c in complexes for (p1, p2) in it.combinations(complexes[c], 2)}

    return complexes_pairs


def get_complexes_dict(organism='Human', corum_file='./files/allComplexesCore.csv'):
    complexes = read_csv(corum_file, sep=';').dropna(subset=['Complex name'])

    complexes = complexes[complexes['organism'] == organism]

    complexes = {c: {x.replace('(', '').replace(')', '') for p in complexes.loc[complexes['Complex id'] == c, 'subunits (UniProt IDs)'] for x in p.split(',')} for c in complexes['Complex id']}

    return complexes


def get_complexes_name(organism='Human', corum_file='./files/allComplexesCore.csv'):
    complexes = read_csv(corum_file, sep=';').dropna(subset=['Complex name'])

    complexes = complexes[complexes['organism'] == organism]

    complexes_name = complexes.groupby('Complex id')['Complex name'].first().to_dict()

    return complexes_name
