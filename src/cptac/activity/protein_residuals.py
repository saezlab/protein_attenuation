import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from pandas import read_csv, DataFrame, concat, Series, pivot_table
from cptac import wd


# -- Correlated protein pairs
p_pairs = read_csv('%s/tables/top_correlated_protein_pairs.csv' % wd)
p_pairs = p_pairs[p_pairs['score'] > .5]
p_pairs = {(p1, p2) for p1, p2 in p_pairs[['p1', 'p2']].values}
print len(p_pairs)


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


# -- Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print proteomics


# -- Estimate residuals
# p = 'LAMB1'
p_residuals, p_predicted, p_lm = {}, {}, {}
for p in network_i.vs['name']:
    if p in proteomics.index:
        c_proteins = set(network_i.vs[network_i.neighborhood(p)]['name']).difference({p}).intersection(proteomics.index)

        if len(c_proteins) > 0:
            y = proteomics.ix[p].dropna()

            x = sm.add_constant(proteomics.ix[c_proteins, y.index].T).dropna()

            if len(x) > 1:
                y = y.ix[x.index]

                lm = sm.OLS(y, x).fit()
                print lm.summary()

                p_residuals[p] = lm.resid
                p_lm[p] = lm

p_residuals = DataFrame(p_residuals).T
p_residuals.to_csv('%s/tables/protein_residuals.csv' % wd)
print p_residuals
