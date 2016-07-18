import igraph
import numpy as np
import seaborn as sns
import itertools as it
import statsmodels.api as sm
import matplotlib.pyplot as plt
from mtkirc.utils import gkn
from pandas.stats.misc import zscore
from cptac.utils import get_cancer_genes
from pymist.utils.stringdb import get_stringdb
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pandas import read_csv, DataFrame, concat, Series, pivot_table
from cptac import wd, palette, palette_cnv, palette_cnv_number, default_color


# -- Protein interactions
uniprot = read_uniprot_genename()

# CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print len(corum)

# String
string = get_stringdb(900)
string = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in string for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print len(string)

# Intersection
p_pairs = corum.intersection(string)
print len(p_pairs)


# -- Proteomics
proteomics = read_csv('%s/data/pancan_preprocessed.csv' % wd, index_col=0)
proteomics.columns = [i[:15] for i in proteomics]

remove_samples = {i for i in set(proteomics) if proteomics.loc[:, [i]].shape[1] == 2 and proteomics.loc[:, [i]].corr().ix[0, 1] < .4}
proteomics = proteomics.drop(remove_samples, axis=1)

proteomics = DataFrame({i: proteomics.loc[:, [i]].mean(1) for i in set(proteomics)})
print proteomics


# -- Create network
network_i = igraph.Graph(directed=False)

# Add nodes
network_i.add_vertices(list({i for p in p_pairs for i in p}))

# Add edges
network_i.add_edges([(p1, p2) for p1, p2 in p_pairs])

# Remove duplicated edges and self-loops
network_i = network_i.simplify(True, True, 'first')
# print DataFrame(Series(dict(zip(*(network_i.vs['name'], network_i.degree())))).sort_values())
print network_i.summary()


# -- Estimate residuals
# p = 'LAMB1'
p_residuals, p_predicted, p_lm = {}, {}, {}
for p in network_i.vs['name']:
    if p in proteomics.index:
        c_proteins = set(network_i.vs[network_i.neighborhood(p)]['name']).difference({p}).intersection(proteomics.index)

        if len(c_proteins) > 0:
            y = proteomics.ix[p].dropna()

            x = sm.add_constant(proteomics.ix[c_proteins, y.index].T).dropna()

            if len(x) > 0:
                y = y.ix[x.index]

                lm = sm.OLS(y, x).fit()
                print lm.summary()

                p_residuals[p] = lm.resid
                p_lm[p] = lm

p_residuals = DataFrame(p_residuals).T
p_residuals.to_csv('%s/tables/protein_residuals.csv' % wd)
print p_residuals
