#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pydot
import igraph
import pickle
import itertools as it
from pandas import DataFrame, read_csv
from protein_attenuation.utils import jaccard, read_uniprot_genename, get_complexes_name, get_complexes_dict


# -- Import proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
proteomics_dict = {s: proteomics[s].dropna().to_dict() for s in proteomics}

samples = set(proteomics)
proteins = set(proteomics.index)


# -- Protein complexes
uniprot, corum_n = read_uniprot_genename(), get_complexes_name()

# dict
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if 1 < len(corum[k])}


# -- Simplify corum complexes sets
corum_jacc = DataFrame([{'c1': c1, 'c2': c2, 'j': jaccard(corum[c1], corum[c2])} for c1, c2 in it.product(corum.keys(), corum.keys())])

# Build network
network_i = igraph.Graph(directed=False)

# Network lists
edges = [(str(px), str(py)) for px, py in corum_jacc.loc[corum_jacc['j'] >= 1., ['c1', 'c2']].values]
vertices = list({p for px, py in edges for p in (px, py)})

# Add nodes
network_i.add_vertices(vertices)

# Add edges
network_i.add_edges(edges)

# Simplify
network_i.simplify(loops=False)

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

with open('./tables/corum_dict_non_redundant.pickle', 'wb') as handle:
    pickle.dump(corum_s, handle, protocol=pickle.HIGHEST_PROTOCOL)
print '[INFO] CORUM complexes redundancy removed (python dictionary): ' + './tables/corum_dict_non_redundant.pickle'
