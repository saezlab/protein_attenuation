#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pydot
import igraph
import numpy as np
import seaborn as sns
import itertools as it
import matplotlib.pyplot as plt
import statsmodels.api as sm
from cptac.utils import jaccard
from cptac import palette_cnv_number
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index
from lifelines.statistics import logrank_test
from statsmodels.stats.weightstats import CompareMeans, DescrStatsW
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LinearRegression
from statsmodels.duration.hazard_regression import PHReg
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict, get_complexes_pairs
from pandas import DataFrame, Series, read_csv, pivot_table, concat


# -- Imports
# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
proteomics_dict = {s: proteomics[s].dropna().to_dict() for s in proteomics}
print 'proteomics', proteomics.shape

# Clinical data
clinical = read_csv('./data/tcga_clinical.csv', index_col=0).dropna(subset=['time', 'status'])
clinical = clinical[(clinical['admin.disease_code'] == 'ov') & (clinical['time'] < (365 * 10))]
print 'clinical', clinical.shape


# -- Overlap
samples, proteins = set(clinical.index).intersection(proteomics), set(proteomics.index)
print 'samples', 'proteins', len(samples), len(proteins)


# -- Protein complexes
uniprot = read_uniprot_genename()

corum_n = get_complexes_name()

# dict
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if 1 < len(corum[k])}
print 'corum', len(corum)

# pairs
corum_pairs = get_complexes_pairs()
corum_pairs = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum_pairs for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print 'corum', len(corum_pairs)


# -- Tumour suppressors and oncogenes
interactions = {(px, py) for px, py in read_csv('./tables/network-112015_symbol.txt', sep='\t').dropna()[['V4', 'V5']].values}
print 'interactions', len(interactions)

cgenes = read_csv('./tables/Cancer5000.oncodriveROLE.0.3-0.7.txt', sep='\t').append(read_csv('./tables/HCD.oncodriveROLE.0.3-0.7.txt', sep='\t'))
cgenes = {
    'suppressor': set(cgenes[cgenes['oncodriveROLE'] == 'Loss of function']['SYM']),
    'oncogene': set(cgenes[cgenes['oncodriveROLE'] == 'Activating']['SYM'])
}

cgenes_omnipath = {
    'suppressor': {(px, py) for px, py in interactions if px in cgenes['suppressor']},
    'oncogene': {(px, py) for px, py in interactions if px in cgenes['oncogene']}
}
print 'cgenes_omnipath', len(cgenes_omnipath['suppressor']), len(cgenes_omnipath['oncogene'])

cgenes_corum = {
    'suppressor': {(px, py) for px, py in corum_pairs if px in cgenes['suppressor'] and (px, py) in cgenes_omnipath['suppressor']},
    'oncogene': {(px, py) for px, py in corum_pairs if px in cgenes['oncogene'] and (px, py) in cgenes_omnipath['suppressor']}
}
print 'cgenes_corum', len(cgenes_corum['suppressor']), len(cgenes_corum['oncogene'])


# -- Simplify corum complexes sets
corum_jacc = DataFrame([{'c1': c1, 'c2': c2, 'j': jaccard(corum[c1], corum[c2])} for c1, c2 in it.product(corum.keys(), corum.keys())])
print corum_jacc.sort('j').tail()

# Build network
network_i = igraph.Graph(directed=False)

# Network lists
edges = [(str(px), str(py)) for px, py in corum_jacc.loc[corum_jacc['j'] > .9, ['c1', 'c2']].values]
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

graph.set_graph_defaults(packMode='clust', pack='true')
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

graph.write_pdf('./reports/corum_jaccard_network.pdf')

# Connected components
components = [network_i.vs[c]['name'] for c in network_i.components()]

# Simplify corum
corum_s = {':'.join(c): {p for i in c for p in corum[int(i)]} for c in components}
print 'corum_s', len(corum_s)


# -- Estimate complex activity
# ts_complex = {px for px, py in cgenes_corum['suppressor']}
ts_complex = {c: corum_s[c] for c in corum_s if len(corum_s[c].intersection(cgenes['suppressor'])) > 0}


def ztest_complex(s, c, df):
    x1 = [df[s][k] for k in df[s] if k in ts_complex[c]]
    x2 = [df[s][k] for k in df[s] if k not in ts_complex[c]]

    if len(x1) > 1:
        stat = CompareMeans(DescrStatsW(x1), DescrStatsW(x2))

        z_larger, p_larger = stat.ztest_ind(alternative='larger')
        z_smaller, p_smaller = stat.ztest_ind(alternative='smaller')

        z, p = z_larger, p_larger if p_larger < p_smaller else p_smaller

        res = {
            'sample': s, 'complex': c, 'name': corum_n[int(c)] if ':' not in c else c,
            'z': z, 'pval': p, 'mean': np.mean(x1), 'targets': len(x1)
        }

        return res

c_activity = [ztest_complex(s, c, proteomics_dict) for s in samples for c in ts_complex]
c_activity = DataFrame([i for i in c_activity if i])
c_activity['fdr'] = multipletests(c_activity['pval'],  method='fdr_bh')[1]
print c_activity.sort('fdr')

# Build matrix
c_activity_matrix = pivot_table(c_activity, index='complex', columns='sample', values='z', fill_value=np.nan)
print 'c_activity_matrix', c_activity_matrix.shape


# # -- Regulatory interactions
# ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
# ppairs_trans = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
# print len(ppairs_trans)
#
# ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
# ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
# ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans for px, py in ppairs_cnv[['px', 'py']].values]]
# print ppairs_cnv.sort('fdr')
#
# associations = {(px, py) for px, py in ppairs_cnv[['px', 'py']].values}
# print len(associations)


# -- Logrank test
# px = '2470'
surv = {}
for px in c_activity_matrix.index:
    df = c_activity_matrix.ix[px, samples].dropna()

    samples_up = set(df[df > (df.mean() + df.std())].index)
    samples_dw = set(df[df < (df.mean() - df.std())].index)
    samples_bg = set(df.index).difference(samples_up.union(samples_dw))

    # samples_up = set(df[df > 0].index)
    # samples_dw = set(df[df < 0].index)
    # samples_bg = set(samples).difference(samples_up.union(samples_dw))

    # if len(samples_up) >= len(samples) * .1 and len(samples_dw) >= len(samples) * .1:
    logrank = logrank_test(
        clinical.ix[samples_up, 'time'], clinical.ix[samples_dw, 'time'],
        clinical.ix[samples_up, 'status'], clinical.ix[samples_dw, 'status']
    )
    print logrank

    res = {
        'pval': logrank.p_value, 't': logrank.test_statistic, 'name': corum_n[int(px)] if ':' not in px else px, 'len': len(corum_s[px].intersection(proteins)),
        'len_up': len(samples_up), 'len_down': len(samples_dw)
    }

    surv[px] = res

surv = DataFrame(surv).T
surv['fdr'] = multipletests(surv['pval'], method='fdr_bh')[1]
print surv.sort(['fdr', 'pval'])

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
ax = plt.subplot(111)

kmf = KaplanMeierFitter()

kmf.fit(clinical.ix[samples_up, 'time'], event_observed=clinical.ix[samples_up, 'status'], label='Up (N=%d)' % len(samples_up))
kmf.plot(ci_force_lines=False, color=palette_cnv_number[1], ax=ax, show_censors=True, censor_styles={'ms': 5, 'marker': 'x'}, lw=1., ci_show=False)

kmf.fit(clinical.ix[samples_bg, 'time'], event_observed=clinical.ix[samples_bg, 'status'], label='Bkg (N=%d)' % len(samples_bg))
kmf.plot(ci_force_lines=False, color=palette_cnv_number[0], ax=ax, show_censors=True, censor_styles={'ms': 5, 'marker': 'x'}, lw=1., ci_show=False)

kmf.fit(clinical.ix[samples_dw, 'time'], event_observed=clinical.ix[samples_dw, 'status'], label='Down (N=%d)' % len(samples_dw))
kmf.plot(ci_force_lines=False, color=palette_cnv_number[-1], ax=ax, show_censors=True, censor_styles={'ms': 5, 'marker': 'x'}, lw=1., ci_show=False)

sns.despine(ax=ax)

ax.set_title('%s\nLogrank test: %.2e' % (corum_n[int(px)] if ':' not in px else px, logrank.p_value))
ax.set_xlabel('Timeline (days)')
ax.set_ylabel('Survival fraction')

plt.gcf().set_size_inches(4, 2)
plt.savefig('./reports/survival_protein_regulators.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
