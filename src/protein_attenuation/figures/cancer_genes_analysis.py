#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pydot
import igraph
import numpy as np
import itertools as it
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from protein_attenuation.utils import jaccard
from scipy.stats.stats import pearsonr, ttest_ind
from matplotlib.gridspec import GridSpec
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.weightstats import CompareMeans, DescrStatsW
from protein_attenuation import palette, default_color, palette_cnv_number
from protein_attenuation.utils import log_likelihood, f_statistic, r_squared
from sklearn.linear_model import LinearRegression
from matplotlib_venn import venn3, venn3_circles
from statsmodels.stats.multitest import multipletests
from lifelines.statistics import logrank_test
from pandas import DataFrame, Series, read_csv, concat, pivot_table
from pymist.utils.corumdb import get_complexes_pairs, get_complexes_name, get_complexes_dict
from pymist.utils.stringdb import get_stringdb
from pymist.utils.map_peptide_sequence import read_uniprot_genename, read_fasta


# -- Imports
# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)
print 'cnv', cnv.shape

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
proteomics_dict = {s: proteomics[s].dropna().to_dict() for s in proteomics}
print 'proteomics', proteomics.shape

# Clinical data
clinical = read_csv('./data/tcga_clinical.csv', index_col=0).dropna(subset=['time', 'status'])
clinical = clinical[(clinical['admin.disease_code'] == 'ov') & (clinical['time'] < (365 * 10))]
print 'clinical', clinical.shape


# -- Overlap
samples = set(cnv).intersection(proteomics).intersection(transcriptomics).intersection(clinical.index)
proteins = set(transcriptomics.index).intersection(proteomics.index)
print 'samples', len(samples)


# -- Protein complexes
uniprot = read_uniprot_genename()

corum_n = get_complexes_name()

# dict
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if 1 < len(corum[k])}
print 'corum', len(corum)


# -- Regulatory interactions
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')

associations = {(px, py) for px, py in ppairs_cnv[['px', 'py']].values}
print len(associations)


# -- Tumour suppressors and oncogenes
cgenes = read_csv('./tables/Cancer5000.oncodriveROLE.0.3-0.7.txt', sep='\t').append(read_csv('./tables/HCD.oncodriveROLE.0.3-0.7.txt', sep='\t'))
cgenes = {
    'suppressor': set(cgenes[cgenes['oncodriveROLE'] == 'Loss of function']['SYM']),
    'oncogene': set(cgenes[cgenes['oncodriveROLE'] == 'Activating']['SYM'])
}
print 'cgenes', len(cgenes)


ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
print ppairs_trans[[px in cgenes['suppressor'] for px, fdr in ppairs_trans[['px', 'fdr']].values]].sort(['fdr', 'f_pval'])

# --
# c = 2721
res = {}
for c in corum:
    if len(corum[c].intersection(cgenes['suppressor'])) > 0 and len(corum[c].intersection(set(ppairs_cnv['px']))) > 0:

        x = cnv.ix[corum[c].intersection(cgenes['suppressor']), samples].T
        y = proteomics.ix[corum[c].intersection(cgenes['suppressor']), samples].T


        x_ts = (x_ts <= -1).any(1).astype(int)

        x_ts_px = cnv.ix[corum[c].intersection(cgenes['suppressor']).union(corum[c].intersection(set(ppairs_cnv['px']))), samples].T
        x_ts_px = (x_ts_px <= -1).any(1).astype(int)

        if len(x_ts[x_ts == 1]) > 2 and len(x_ts_px[x_ts_px == 1]) > 2:
            logrank_ts = logrank_test(
                clinical.ix[x_ts[x_ts == 1].index, 'time'], clinical.ix[x_ts[x_ts != 1].index, 'time'],
                clinical.ix[x_ts[x_ts == 1].index, 'status'], clinical.ix[x_ts[x_ts != 1].index, 'status']
            )

            logrank_ts_px = logrank_test(
                clinical.ix[x_ts_px[x_ts_px == 1].index, 'time'], clinical.ix[x_ts_px[x_ts_px != 1].index, 'time'],
                clinical.ix[x_ts_px[x_ts_px == 1].index, 'status'], clinical.ix[x_ts_px[x_ts_px != 1].index, 'status']
            )

            res[c] = {'ts': logrank_ts.p_value, 'ts_px': logrank_ts_px.p_value}

res = DataFrame(res).T
print res

plt.scatter(res['ts'], res['ts_px'])


# # -- Protein complexes
# uniprot = read_uniprot_genename()
#
# corum_n = get_complexes_name()
#
# # dict
# corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
# corum = {k: corum[k] for k in corum if 1 < len(corum[k])}
# print 'corum', len(corum)


# # -- Simplify corum complexes sets
# corum_jacc = DataFrame([{'c1': c1, 'c2': c2, 'j': jaccard(corum[c1], corum[c2])} for c1, c2 in it.product(corum.keys(), corum.keys())])
# print corum_jacc.sort('j').tail()
#
# # Build network
# network_i = igraph.Graph(directed=False)
#
# # Network lists
# edges = [(str(px), str(py)) for px, py in corum_jacc.loc[corum_jacc['j'] > .9, ['c1', 'c2']].values]
# vertices = list({p for px, py in edges for p in (px, py)})
#
# # Add nodes
# network_i.add_vertices(vertices)
# print network_i.summary()
#
# # Add edges
# network_i.add_edges(edges)
# print network_i.summary()
#
# # Simplify
# network_i.simplify(loops=False)
# print network_i.summary()
#
# # Print
# graph = pydot.Dot(graph_type='graph')
#
# graph.set_graph_defaults(packMode='clust', pack='true')
# graph.set_node_defaults(fontcolor='white', penwidth='5', fillcolor='#CCCCCC', width='1', height='1', fontsize='20', fontname='sans-serif')
# graph.set_edge_defaults(color='#CCCCCC', arrowhead='vee', penwidth='2.')
#
# for e in network_i.es:
#     if e.source != e.target:
#         source = pydot.Node(network_i.vs[e.source]['name'], style='filled', shape='ellipse', penwidth='0')
#         target = pydot.Node(network_i.vs[e.target]['name'], style='filled', shape='ellipse', penwidth='0')
#
#         graph.add_node(source)
#         graph.add_node(target)
#
#         edge = pydot.Edge(source, target)
#         graph.add_edge(edge)
#
# graph.write_pdf('./reports/corum_jaccard_network.pdf')
#
# # Connected components
# components = [network_i.vs[c]['name'] for c in network_i.components()]
#
# # Simplify corum
# corum_s = {':'.join(c): {p for i in c for p in corum[int(i)]} for c in components}
# print 'corum_s', len(corum_s)
#
#
# # -- Estimate complex activity
# # ts_complex = {px for px, py in cgenes_corum['suppressor']}
# ts_complex = {c: corum_s[c] for c in corum_s if len(corum_s[c].intersection(cgenes['suppressor'])) > 0}
#
#
# def ztest_complex(s, c, df):
#     x1 = [df[s][k] for k in df[s] if k in ts_complex[c]]
#     x2 = [df[s][k] for k in df[s] if k not in ts_complex[c]]
#
#     if len(x1) > 1:
#         stat = CompareMeans(DescrStatsW(x1), DescrStatsW(x2))
#
#         z_larger, p_larger = stat.ztest_ind(alternative='larger')
#         z_smaller, p_smaller = stat.ztest_ind(alternative='smaller')
#
#         z, p = z_larger, p_larger if p_larger < p_smaller else p_smaller
#
#         res = {
#             'sample': s, 'complex': c, 'name': corum_n[int(c)] if ':' not in c else c,
#             'z': z, 'pval': p, 'mean': np.mean(x1), 'targets': len(x1)
#         }
#
#         return res
#
# c_activity = [ztest_complex(s, c, proteomics_dict) for s in samples for c in ts_complex]
# c_activity = DataFrame([i for i in c_activity if i])
# c_activity['fdr'] = multipletests(c_activity['pval'],  method='fdr_bh')[1]
# print c_activity.sort('fdr')
#
# # Build matrix
# c_activity_matrix = pivot_table(c_activity, index='complex', columns='sample', values='z', fill_value=np.nan)
# print 'c_activity_matrix', c_activity_matrix.shape
#
#
# # -- Protein correlations
# p_correlations = read_csv('./tables/proteins_correlations.csv', index_col=0)
# print 'p_correlations', len(p_correlations)
#
#
# # -- Plot scatter of correlations
# ax_min = np.min([p_correlations['CNV_Transcriptomics'].min() * 1.10, p_correlations['Transcriptomics_Proteomics'].min() * 1.10])
# ax_max = np.min([p_correlations['CNV_Transcriptomics'].max() * 1.10, p_correlations['Transcriptomics_Proteomics'].max() * 1.10])
#
# sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'lines.linewidth': .75})
# g = sns.jointplot(
#     'CNV_Transcriptomics', 'CNV_Proteomics', p_correlations, 'scatter', color='#808080', xlim=[ax_min, ax_max], ylim=[ax_min, ax_max],
#     space=0, s=15, edgecolor='w', linewidth=.1, marginal_kws={'hist': False, 'rug': False}, stat_func=None, alpha=.3
# )
# g.plot_marginals(sns.kdeplot, shade=True, color='#595959', lw=.3)
#
# g.ax_joint.axhline(0, ls='-', lw=0.1, c='black', alpha=.3)
# g.ax_joint.axvline(0, ls='-', lw=0.1, c='black', alpha=.3)
# g.ax_joint.plot([ax_min, ax_max], [ax_min, ax_max], 'k--', lw=.3)
#
# g.x = p_correlations['CNV_Transcriptomics']
# g.y = p_correlations['CNV_Proteomics']
# g.plot_joint(sns.kdeplot, cmap=sns.light_palette('#595959', as_cmap=True), legend=False, shade=False, shade_lowest=False, n_levels=9, alpha=.8, lw=.1)
#
# g.x = p_correlations[[i in cgenes['suppressor'] for i in p_correlations.index]]['CNV_Transcriptomics']
# g.y = p_correlations[[i in cgenes['suppressor'] for i in p_correlations.index]]['CNV_Proteomics']
# g.plot_joint(sns.regplot, fit_reg=False, color=palette_cnv_number[-2])
# g.plot_marginals(sns.kdeplot, color=palette_cnv_number[-2], shade=True, legend=False)
#
# g.x = p_correlations[[i in cgenes['oncogene'] for i in p_correlations.index]]['CNV_Transcriptomics']
# g.y = p_correlations[[i in cgenes['oncogene'] for i in p_correlations.index]]['CNV_Proteomics']
# g.plot_joint(sns.regplot, fit_reg=False, color=palette_cnv_number[2])
# g.plot_marginals(sns.kdeplot, color=palette_cnv_number[2], shade=True, legend=False)
#
# plt.gcf().set_size_inches(3, 3)
#
# g.set_axis_labels('CNV ~ Transcriptomics', 'CNV ~ Proteomics')
# plt.savefig('./reports/correlation_difference_lmplot_corr_cancer_genes.pdf', bbox_inches='tight')
# plt.savefig('./reports/correlation_difference_lmplot_corr_cancer_genes.png', bbox_inches='tight', dpi=300)
# plt.close('all')
# print '[INFO] Plot done'
#
#
# # -- Protein complexes interactions
# uniprot = read_uniprot_genename()
#
# # string = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in get_stringdb(900) for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
# # print 'string', len(string)
#
# corum = get_complexes_pairs()
# corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
# print 'corum', len(corum)
#
#
# # -- Overlap
# samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
# proteins = {p for (px, py) in corum for p in [px, py]}.intersection(proteomics.index)
# print 'samples', len(samples)
#
#
# # --
# # g = 'ZW10'
# def protein_residual(g):
#     y = proteomics.ix[g, samples].dropna()
#     x = transcriptomics.ix[[g], y.index].T
#
#     lm = LinearRegression().fit(x, y)
#     y_ = y - lm.coef_[0] * x[g] - lm.intercept_
#
#     return y_
# residuals = DataFrame({g: protein_residual(g) for g in set(transcriptomics.index).intersection(proteomics.index)}).T
# print 'residuals', residuals.shape
#
#
# # -- Tumour suppressors and oncogenes
# cgenes = read_csv('./tables/Cancer5000.oncodriveROLE.0.3-0.7.txt', sep='\t').append(read_csv('./tables/HCD.oncodriveROLE.0.3-0.7.txt', sep='\t'))
# cgenes = {
#     'suppressor': set(cgenes[cgenes['oncodriveROLE'] == 'Loss of function']['SYM']),
#     'oncogene': set(cgenes[cgenes['oncodriveROLE'] == 'Activating']['SYM'])
# }
# print 'cgenes', len(cgenes)
#
# # interactions = {(px, py) for px, py in read_csv('./tables/network-112015_symbol.txt', sep='\t').dropna()[['V4', 'V5']].values if px in cgenes['suppressor'] or px in cgenes['oncogene']}
# interactions = {(px, py) for px, py in read_csv('./tables/network-112015_symbol.txt', sep='\t').dropna()[['V4', 'V5']].values if px in cgenes['suppressor']}
# print 'interactions', len(interactions)
#
# sources, targets = {px for px, py in interactions}, {py for px, py in interactions}
# print 'sources', 'targets', len(sources), len(targets)
#
#
# # -- Regressions: Py Residuals ~ Px CNV
# # px, py = 'SMARCA4', 'ZEB1'
# def regressions(px, py):
#     if px in cnv.index and py in residuals.index and px != py:
#         # Protein measurements
#         y = residuals.ix[py, samples].dropna()
#         x = cnv.ix[[px], y.index].T.dropna()
#
#         y = y.ix[x.index]
#
#         # Fit models
#         lm = LinearRegression().fit(x, y)
#
#         # Predict
#         y_true, y_pred = y.copy(), Series(dict(zip(*(x.index, lm.predict(x)))))
#
#         # Log likelihood
#         l_lm = log_likelihood(y_true, y_pred)
#
#         # F-statistic
#         f, f_pval = f_statistic(y_true, y_pred, len(y), x.shape[1])
#
#         # R-squared
#         r = r_squared(y_true, y_pred)
#
#         res = {
#             'px': px, 'py': py, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta': lm.coef_[0]
#         }
#
#         print 'Px (%s), Py (%s): Rsquared: %.2f, F: %.2f, F pval: %.2e' % (px, py, res['rsquared'], res['f'], res['f_pval'])
#         # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()
#
#         return res
#
# ppairs = [regressions(px, py) for px in sources for py in targets]
# ppairs = DataFrame([i for i in ppairs if i])
# ppairs['fdr'] = multipletests(ppairs['f_pval'], method='fdr_bh')[1]
# ppairs['type'] = ['suppressor' if i in cgenes['suppressor'] else 'oncogene' for i in ppairs['px']]
# print ppairs.sort(['fdr', 'f_pval'])
#
#
# ppairs[ppairs['fdr'] < .05].groupby('px')['py'].agg(lambda x: set(x)).to_dict()
#
# sns.boxplot(x='type', y='beta', data=ppairs)
# sns.boxplot(x='type', y='beta', data=ppairs[ppairs['fdr'] < .05])
