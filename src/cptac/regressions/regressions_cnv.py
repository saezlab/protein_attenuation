import pydot
import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from scipy.stats.stats import pearsonr, ttest_ind
from matplotlib.gridspec import GridSpec
from statsmodels.stats.weightstats import ztest
from cptac import wd, palette, default_color, palette_cnv_number
from cptac.utils import log_likelihood, f_statistic, r_squared
from sklearn.linear_model import LinearRegression
from matplotlib_venn import venn3, venn3_circles
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename, read_fasta


# -- Imports
# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
print 'cnv', cnv.shape

# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print 'proteomics', proteomics.shape

# -- Overlap
genes = set(cnv.index).intersection(proteomics.index).intersection(transcriptomics.index)
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
print 'genes', 'samples', len(genes), len(samples)


# -- Protein complexes interactions
uniprot = read_uniprot_genename()
uniprot_fasta = read_fasta()

corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
corum = {(p1, p2) for p1, p2 in corum if p1 in genes and p2 in genes}
print 'corum', len(corum)


# --
# g = 'ZW10'
def protein_residual(g):
    y = proteomics.ix[g, samples].dropna()
    x = transcriptomics.ix[[g], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[g] - lm.intercept_

    return y_
residuals = DataFrame({g: protein_residual(g) for g in genes}).T
print 'residuals', residuals.shape


# -- Regressions: Py Residuals ~ Px CNV
# px, py = 'ARID1A', 'DPF2'
def regressions(px, py):
    # Protein measurements
    y = residuals.ix[py, samples].dropna()
    x = cnv.ix[[px], y.index].T

    # Fit models
    lm = LinearRegression().fit(x, y)

    # Predict
    y_true, y_pred = y.copy(), Series(dict(zip(*(x.index, lm.predict(x)))))

    # Log likelihood
    l_lm = log_likelihood(y_true, y_pred)

    # F-statistic
    f, f_pval = f_statistic(y_true, y_pred, len(y), x.shape[1])

    # R-squared
    r = r_squared(y_true, y_pred)

    res = {
        'px': px, 'py': py, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm
    }

    print 'Px (%s), Py (%s): Rsquared: %.2f, F: %.2f, F pval: %.2e' % (px, py, res['rsquared'], res['f'], res['f_pval'])
    # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()

    return res

ppairs = DataFrame([regressions(px, py) for px, py in corum])
ppairs['fdr'] = multipletests(ppairs['f_pval'], method='fdr_bh')[1]
ppairs.to_csv('%s/tables/ppairs_cnv_regulation_all.csv' % wd, index=False)
# ppairs = read_csv('%s/tables/ppairs_cnv_regulation_all.csv' % wd)
print ppairs.sort('fdr')


# -- QQ-plot
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})

plot_df = DataFrame({
    'x': sorted([-np.log10(np.float(i) / len(ppairs)) for i in np.arange(1, len(ppairs) + 1)]),
    'y': sorted(-np.log10(ppairs['f_pval']))
})

g = sns.regplot('x', 'y', plot_df, fit_reg=False, ci=False, color=default_color, line_kws={'lw': .3})
g.set_xlim(0)
g.set_ylim(0)
plt.plot(plt.xlim(), plt.xlim(), 'k--', lw=.3)
plt.xlabel('Theoretical -log(P)')
plt.ylabel('Observed -log(P)')
plt.title('Residuals ~ Copy Number')
sns.despine(trim=True)
plt.gcf().set_size_inches(3, 3)
plt.savefig('%s/reports/ppairs_cnv_regulation_qqplot.png' % wd, bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Plot done'


# -- Create network
network_i = igraph.Graph(directed=True)

# Initialise network lists
edges = [(px, py) for px, py in ppairs[ppairs['fdr'] < .05][['px', 'py']].values]
vertices = list({p for px, py in edges for p in (px, py)})

# Add nodes
network_i.add_vertices(vertices)
print network_i.summary()

# Add edges
network_i.add_edges(edges)
print network_i.summary()

# Draw network
graph = pydot.Dot(graph_type='digraph', rankdir='LR')

graph.set_node_defaults(fontcolor='white', penwidth='3', fillcolor='#CCCCCC')
graph.set_edge_defaults(color='#CCCCCC', arrowhead='vee')

for edge in network_i.es:
    source_id, target_id = network_i.vs[[edge.source, edge.target]]['name']

    source = pydot.Node(source_id, style='filled', shape='ellipse', penwidth='0')
    target = pydot.Node(target_id, style='filled', shape='ellipse', penwidth='0')

    graph.add_node(source)
    graph.add_node(target)

    edge = pydot.Edge(source, target)
    graph.add_edge(edge)

graph.write_pdf('%s/reports/ppairs_cnv_regulation_network.pdf' % wd)
print '[INFO] Network PDF exported'


# -- Scatter
def ppair_correlation(px, py):
    x, y = zip(*proteomics.ix[[px, py]].T.dropna().values)
    return pearsonr(x, y)

ppairs_signif = ppairs[ppairs['fdr'] < .05].sort('fdr')
ppairs_signif['cor'] = [ppair_correlation(px, py)[0] for px, py in ppairs_signif[['px', 'py']].values]
ppairs_signif.to_csv('%s/tables/ppairs_cnv_regulation.csv' % wd, index=False)

# px, py = 'COG3', 'COG2'
gs, pos = GridSpec(len(ppairs_signif), 2, hspace=.5), 0
for px, py in ppairs_signif[['px', 'py']].values:
    #
    ax = plt.subplot(gs[pos])

    y = residuals.ix[py, samples].dropna()
    x = cnv.ix[px, y.index]
    plot_df = concat([x, y], axis=1).dropna()

    sns.regplot(plot_df[px], plot_df[py], ax=ax, color=default_color, fit_reg=True, scatter=True, truncate=True)
    for c in [0, -1, 1, -2, 2]:
        sns.regplot(plot_df[plot_df[px] == c][px], plot_df[plot_df[px] == c][py], ax=ax, color=palette_cnv_number[c], fit_reg=False, truncate=True)
    sns.despine(ax=ax)
    ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_xlabel('%s (copy number)' % px)
    ax.set_ylabel('%s (residuals)' % py)
    ax.set_title('Pearson\'s r: %.2f, p-value: %.2e' % pearsonr(x, y))
    ax.set_ylim(plot_df[py].min() * 1.05, plot_df[py].max() * 1.05)

    #
    ax = plt.subplot(gs[pos + 1])

    y = proteomics.ix[py, samples].dropna()
    x = proteomics.ix[px, y.index]
    plot_df = concat([x, y, cnv.ix[px, y.index].rename('cnv')], axis=1).dropna()

    sns.regplot(plot_df[px], plot_df[py], ax=ax, color=default_color, fit_reg=True, scatter=True, truncate=True)
    for c in [0, -1, 1, -2, 2]:
        sns.regplot(plot_df[plot_df['cnv'] == c][px], plot_df[plot_df['cnv'] == c][py], ax=ax, color=palette_cnv_number[c], fit_reg=False, truncate=True)
    sns.despine(ax=ax)
    ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_xlabel('%s (proteomics)' % px)
    ax.set_ylabel('%s (proteomics)' % py)
    ax.set_title('Pearson\'s r: %.2f, p-value: %.2e' % pearsonr(plot_df[px], plot_df[py]))
    ax.set_ylim(plot_df[py].min() * 1.05, plot_df[py].max() * 1.05)

    pos += 2

plt.gcf().set_size_inches(4, 2 * len(ppairs_signif))
plt.savefig('%s/reports/ppairs_cnv_regulation_scatter.png' % wd, bbox_inches='tight', dpi=150)
plt.close('all')
print '[INFO] Plot done'


# -- Protein sequence length
p_info = read_csv('%s/files/expasy_mol_weight.txt' % wd, sep='\t')
p_info = p_info[p_info['pI'] != 'UNDEFINED']
p_info = p_info.groupby('id').max()
p_info['length'] = [len(uniprot_fasta[i]) for i in p_info.index]
p_info.index = [uniprot[i][0] for i in p_info.index]

plot_df = DataFrame([{'type': t, 'protein': p, 'signif': int(fdr < .05), 'feature': f, 'value': p_info.ix[p, f]} for px, py, fdr in ppairs[['px', 'py', 'fdr']].values for t, p in [('px', px), ('py', py)] for f in ['length'] if p in p_info.index])
plot_df = plot_df[plot_df['signif'] == 1]
print plot_df.sort('signif')

ttest, pval = ttest_ind(plot_df[plot_df['type'] == 'px']['value'], plot_df[plot_df['type'] == 'py']['value'])
print 'ttest', 'pval', ttest, pval

# Boxplot
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
sns.violinplot(x='type', y='value', data=plot_df, linewidth=.3, cut=0, inner='quartile', split=False, color=palette_cnv_number[0])
sns.stripplot(x='type', y='value', data=plot_df, linewidth=.3, jitter=True, edgecolor='white', split=False, color=default_color)
plt.ylim(0)
sns.despine(trim=True)
plt.ylabel('Protein sequence length (number of AA)')
plt.title('Px (CNV) ~ Py (Residuals)\nT-test: %.2f, p-value: %.2e' % (ttest, pval))
plt.gcf().set_size_inches(2, 4)
plt.savefig('%s/reports/protein_pairs_protein_info_boxplots.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'

# Histogram
plot_df = DataFrame([{'px': px, 'py': py, 'len_px': p_info.ix[px, 'length'], 'len_py': p_info.ix[py, 'length'], 'signif': int(fdr < .05)} for px, py, fdr in ppairs[['px', 'py', 'fdr']].values if px in p_info.index and py in p_info.index])
plot_df['diff'] = plot_df['len_px'] - plot_df['len_py']
print plot_df.sort('signif')

sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
for i in [0, 1]:
    values = plot_df[plot_df['signif'] == i]['diff']
    sns.distplot(values, hist=False, kde_kws={'shade': True}, color=palette_cnv_number[i], label='%s (mean = %.2f)' % ('Significant' if i else 'All', np.mean(values)))
plt.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
plt.title('Protein length difference')
plt.xlabel('len(Px) - len(Py)')
sns.despine(trim=True)
plt.legend()
plt.gcf().set_size_inches(4, 2)
plt.savefig('%s/reports/protein_pairs_protein_info_histogram.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'



# # -- Plot
# # Significant associations venn diagram
# sns.set(style='white', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
# fig, gs, pos = plt.figure(figsize=(5, 10)), GridSpec(1, 2, hspace=.3), 0
#
# for c in [-2, 2]:
#     ax = plt.subplot(gs[pos])
#
#     associations = {
#         d: {p for p, t, f in res[d][['protein', 'cnv', 'f_adjpval']].values if c == t and f < .05}
#         for d in res}
#
#     venn3(associations.values(), set_labels=associations.keys(), set_colors=[palette[k] for k in associations])
#     venn3_circles(associations.values(), linestyle='solid', color='white')
#
#     ax.set_title('Depletion' if c == -2 else 'Amplification')
#
#     pos += 1
#
# plt.savefig('%s/reports/regressions_overlap_venn.pdf' % wd, bbox_inches='tight')
# plt.close('all')
# print '[INFO] Done'

