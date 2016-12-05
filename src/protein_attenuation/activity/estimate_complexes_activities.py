#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pickle
import pydot
import igraph
import numpy as np
import itertools as it
from protein_attenuation.utils import jaccard
from statsmodels.stats.weightstats import CompareMeans, DescrStatsW
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
proteomics_dict = {s: proteomics[s].dropna().to_dict() for s in proteomics}
print 'proteomics', proteomics.shape


# -- Overlap
samples = set(proteomics)
proteins = set(proteomics.index)
print 'samples', len(samples)


# -- Protein complexes
uniprot, corum_n = read_uniprot_genename(), get_complexes_name()

# dict
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if 1 < len(corum[k])}
print 'corum', len(corum)


# -- Simplify corum complexes sets
corum_jacc = DataFrame([{'c1': c1, 'c2': c2, 'j': jaccard(corum[c1], corum[c2])} for c1, c2 in it.product(corum.keys(), corum.keys())])
print corum_jacc.sort('j').tail()

# Build network
network_i = igraph.Graph(directed=False)

# Network lists
edges = [(str(px), str(py)) for px, py in corum_jacc.loc[corum_jacc['j'] >= 1., ['c1', 'c2']].values]
vertices = list({p for px, py in edges for p in (px, py)})

# Add nodes
network_i.add_vertices(vertices)
print network_i.summary()

# Add edges
network_i.add_edges(edges)
print network_i.summary()

# Simplify
network_i.simplify(loops=False)
print network_i.summary()

# Print
graph = pydot.Dot(graph_type='graph')

graph.set_graph_defaults(packMode='clust', pack='true', dpi='300')
graph.set_node_defaults(fontcolor='white', penwidth='5', fillcolor='#CCCCCC', width='1', height='1', fontsize='20', fontname='sans-serif')
graph.set_edge_defaults(color='#CCCCCC', arrowhead='vee', penwidth='2.')

for e in network_i.es:
    if e.source != e.target:
        source = pydot.Node(network_i.vs[e.source]['name'], style='filled', shape='ellipse', penwidth='0')
        target = pydot.Node(network_i.vs[e.target]['name'], style='filled', shape='ellipse', penwidth='0')

        graph.add_node(source)
        graph.add_node(target)

        edge = pydot.Edge(source, target)
        graph.add_edge(edge)

graph.write_png('./reports/corum_jaccard_network.png')
graph.write_pdf('./reports/corum_jaccard_network.pdf')

# Connected components
components = [network_i.vs[c]['name'] for c in network_i.components()]

# Simplify corum
corum_s = {':'.join(c): {p for i in c for p in corum[int(i)]} for c in components}
print 'corum_s', len(corum_s)

with open('./tables/corum_dict_non_redundant.pickle', 'wb') as handle:
    pickle.dump(corum_s, handle, protocol=pickle.HIGHEST_PROTOCOL)


# -- Estimate complex activity
c_activity = []
for s in samples:
    for c in corum_s:
        x1 = [proteomics_dict[s][k] for k in proteomics_dict[s] if k in corum_s[c]]
        x2 = [proteomics_dict[s][k] for k in proteomics_dict[s] if k not in corum_s[c]]

        if len(x1) > 1:
            stat = CompareMeans(DescrStatsW(x1), DescrStatsW(x2))

            z_larger, p_larger = stat.ztest_ind(alternative='larger')
            z_smaller, p_smaller = stat.ztest_ind(alternative='smaller')

            z, p = z_larger, p_larger if p_larger < p_smaller else p_smaller

            res = {
                'sample': s, 'complex': c, 'name': corum_n[int(c)] if ':' not in c else c,
                'z': z, 'pval': p, 'mean': np.mean(x1), 'targets': len(x1)
            }

            c_activity.append(res)

c_activity = DataFrame(c_activity)
c_activity['fdr'] = multipletests(c_activity['pval'],  method='fdr_bh')[1]
c_activity.to_csv('./tables/protein_complexes_proteomics_activities.csv', index=False)
print c_activity.sort('fdr')
