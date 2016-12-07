#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

from __future__ import division
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from pandas import DataFrame, Series, read_csv
from statsmodels.stats.multitest import multipletests


def slapenrich(mutation_m, pathways, background):
    # Initial variables
    number_bkg_genes, number_samples = len(background), mutation_m.shape[1]

    # R imports: poibin package
    importr('poibin')

    # Pathway stats
    pathways_mutation_status, pathways_probabilities = {}, {}

    for pathway in pathways:
        # Initialise variables
        pathways_mutation_status[pathway], pathways_probabilities[pathway] = {}, {}

        # Overlap pathway genes with mutation matrix genes and background genes
        pathway_genes = pathways[pathway].intersection(mutation_m.index).intersection(background)

        # Pathway specific mutation matrix
        pathway_mutation_m = mutation_m.ix[pathway_genes]

        # Proceed with analysis only if pathway has overlapping genes and at least one gene in the pathway is mutated
        if (len(pathway_genes) > 0) and (pathway_mutation_m.sum().sum() > 0):
            for sample in pathway_mutation_m:
                n_sample_mutations = mutation_m[sample].sum()

                # Check mutational status of the pathway
                pathways_mutation_status[pathway][sample] = int(pathway_mutation_m[sample].sum() > 0)

                # Calculate pathway probablilty of being mutated
                pathways_probabilities[pathway][sample] = 1 - np.exp(-(n_sample_mutations / number_bkg_genes) * len(pathway_genes))

    # Transform pathway status dictionaries into DataFrames
    pathways_mutation_status, pathways_probabilities = DataFrame(pathways_mutation_status), DataFrame(pathways_probabilities)

    # Calculate pathway logOddRatios
    pathway_logoddratios = np.log10(pathways_mutation_status.sum().divide(pathways_probabilities.sum()))

    # Discard pathways with non-overlapping genes
    pathways_ = set(pathway_logoddratios.dropna().index)

    # Calculate pathway p-values
    pathways_pvalues = {}
    for pathway in pathways_:
        # Range between number of samples this pathway was mutated till maximum possible number of mutations
        distribution = range(pathways_mutation_status[pathway].sum(), (number_samples + 1))

        # Pathway probabilities in each sample
        probabilities = list(pathways_probabilities[pathway])

        # Call Poisson Binomial distribution density function
        ro.globalenv['kk'], ro.globalenv['pp'] = distribution, probabilities
        pvalue = sum(ro.r('dpoibin(kk, pp)'))

        # Store result
        pathways_pvalues[pathway] = pvalue
    pathways_pvalues = Series(pathways_pvalues)

    # Adjust p-values
    pathways_fdr = dict(zip(*(pathways_pvalues.index, multipletests(pathways_pvalues, method='fdr_bh')[1])))

    # Export results as dataframe
    results = [(pathway, pathway_logoddratios[pathway], pathways_pvalues[pathway], pathways_fdr[pathway]) for pathway in pathways_]
    results = DataFrame(results, columns=['pathway', 'logoddratios', 'pvalue', 'fdr']).set_index('pathway')

    # Return
    return results, pathways_probabilities


def slapenrich_test(error_threshold=1e-10):
    # Imports
    mutation_m = read_csv('./files/slapenrich_test_genomics.csv', index_col=0)

    background = set(mutation_m.index)

    pathways = read_csv('./files/slapenrich_test_kegg_sets.txt', sep='\t')
    pathways = pathways.groupby('Pathway')['Gene'].agg(lambda x: set(x)).to_dict()

    # Run SlapEnrich python implementation
    results = slapenrich(mutation_m, pathways, background)

    # Import SlapEnrich R results
    results_r = read_csv('./files/slapenrich_test_results_r.csv', index_col=0)

    # Evaluate results
    ratios_error = results_r.ix[results.index, 'logoddratios'].subtract(results['logoddratios']).sum()
    pvalues_error = results_r.ix[results.index, 'pvalue'].subtract(results['pvalue']).sum()

    # Check
    if ratios_error > error_threshold:
        print '[WARNING] Ratios not correct: %.1e' % ratios_error

    elif pvalues_error > error_threshold:
        print '[WARNING] P-values not correct: %.1e' % pvalues_error

    else:
        print '[INFO] Ratios and P-values correctly calculated within %.1e threshold' % error_threshold

    return (ratios_error <= error_threshold) and (pvalues_error <= error_threshold)
