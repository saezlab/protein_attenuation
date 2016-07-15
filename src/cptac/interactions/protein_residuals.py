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


def randomise_matrix(matrix):
    random_df = matrix.copy()
    movers = ~np.isnan(random_df.values)
    random_df.values[movers] = np.random.permutation(random_df.values[movers])
    return random_df

# -- Uniprot
uniprot = read_uniprot_genename()
cancer_genes = get_cancer_genes()

# -- Imports
# CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print len(corum)

# String
string = get_stringdb(900)
string = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in string for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print len(string)

# Intersection
p_pairs = corum.union(string)
print len(p_pairs)

# -- CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
cnv.columns = [i[:15] for i in cnv]

cnv_df = cnv.unstack().reset_index()
cnv_dict = cnv.to_dict()
print cnv

# -- Proteomics
proteomics = read_csv('%s/data/pancan_preprocessed_normalised.csv' % wd, index_col=0)
proteomics = proteomics[proteomics.count(1) > (proteomics.shape[1] * .5)]
annot = read_csv('%s/data/samplesheet.csv' % wd, index_col=0)
print proteomics

# -- Gexp
transcriptomics = read_csv('%s/data/tcga_rnaseq.tsv' % wd, sep='\t', index_col=0)
transcriptomics_annot = Series({i: c for c in transcriptomics for i in annot.index if c[:16] == i[:16]})
transcriptomics = transcriptomics[transcriptomics_annot.values]
transcriptomics.columns = transcriptomics_annot.index
print transcriptomics

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


# -- Import proteomics
brca = read_csv('%s/data/brca_proteomics_processed.csv' % wd, index_col=0)
coread = read_csv('%s/data/coread_proteomics_processed.csv' % wd, index_col=0)
hgsc = read_csv('%s/data/hgsc_proteomics_processed.csv' % wd, index_col=0)

# Concatenate all
pancan = concat([brca, hgsc, coread], axis=1)
pancan = pancan[pancan.count(1) > (pancan.shape[1] * .5)]
print pancan

# -- Covariates
# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()

design = Series(np.concatenate([
    np.repeat('brca', brca.shape[1]),
    np.repeat('hgsc', hgsc.shape[1]),
    np.repeat('coread', coread.shape[1])
]), index=pancan.columns)
design = design.str.get_dummies()

design = concat([design, Series({i: clinical_gender['-'.join(i.split('-')[:4])[:-1].upper()].lower() for i in design.index}).str.get_dummies()], axis=1)
design['age'] = [clinical_age['-'.join(i.split('-')[:4])[:-1].upper()] for i in design.index]
design['tmt'] = np.bitwise_or(design['brca'], design['hgsc'])
design['shotgun'] = design['coread']
design['const'] = 1
print list(design)


# p = 'LAMB1'
p_residuals, p_predicted, p_lm = {}, {}, {}
for p in network_i.vs['name']:
    if p in pancan.index:
    # if p in proteomics.index:
        c_proteins = set(network_i.vs[network_i.neighborhood(p)]['name']).difference({p}).intersection(pancan.index)

        if len(c_proteins) > 0:
            y = pancan.ix[p].dropna()
            # y = proteomics.ix[p].dropna()

            x = concat([pancan.ix[c_proteins].T, design], axis=1).ix[y.index].dropna()
            # x = sm.add_constant(proteomics.ix[c_proteins, y.index].T).dropna()

            if len(x) > 0:
                y = y.ix[x.index]

                lm = sm.OLS(y, x).fit()
                print lm.summary()

                p_residuals[p] = lm.resid
                p_predicted[p] = dict(zip(*(x.index, lm.predict(x))))
                p_lm[p] = lm

p_residuals = DataFrame(p_residuals).T
p_predicted = DataFrame(p_predicted).T
p_residuals.to_csv('%s/tables/protein_residuals.csv' % wd)
print p_residuals

# # --
# cnv_down = {(b, g) for b, g in cnv_df.loc[cnv_df[0] == -2, ['level_0', 'hgnc']].values}
# cnv_up = {(b, g) for b, g in cnv_df.loc[cnv_df[0] == 2, ['level_0', 'hgnc']].values}
#
# plot_df = []
# for t in np.arange(0, 1, .1):
#     t_proteins = {p for p in p_lm if p_lm[p].rsquared_adj > t}
#
#     df_ = DataFrame({i: gkn(p_residuals.ix[i].dropna()).to_dict() for i in t_proteins}).T
#     df_ = df_.unstack().reset_index().dropna()
#     df_.columns = ['barcode', 'protein', 'value']
#
#     cnv_down_mean = np.mean([v for b, p, v in df_[['barcode', 'protein', 'value']].values if (b[:15], p) in cnv_down])
#     cnv_up_mean = np.mean([v for b, p, v in df_[['barcode', 'protein', 'value']].values if (b[:15], p) in cnv_up])
#
#     cnv_diff = cnv_up_mean - cnv_down_mean
#
#     plot_df.append((t, cnv_diff))
#
# plot_df = DataFrame(plot_df, columns=['threshold', 'difference'])
# print plot_df
#
# sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
# g = sns.jointplot(
#     'threshold', 'difference', data=plot_df, color=default_color, annot_kws={'template': 'Pearson\'s r: {val:.2g}, p-value: {p:.1e}', 'loc': 3},
#     joint_kws={'s': 40, 'edgecolor': 'w', 'linewidth': .5}, marginal_kws={'hist': False, 'rug': True}, space=0
# )
# g.set_axis_labels('R-square threshold', 'mean(CNV amplified) - mean(CNV depleted)')
# plt.suptitle('CNV difference by thresholds')
# plt.savefig('%s/reports/pancan_residuals_theshold_cnv_scatter.pdf' % wd, bbox_inches='tight')
# plt.close('all')
# print '[INFO] Done'


# # --
# plot_df = DataFrame([])
# for n, df in [('residuals', p_residuals), ('residuals (random)', randomise_matrix(p_residuals)), ('proteomics', proteomics), ('transcriptomics', transcriptomics)]:
#     df_ = DataFrame({i: gkn(df.ix[i].dropna()).to_dict() for i in df.index}).T
#     df_ = df_.unstack().reset_index().dropna()
#     df_.columns = ['barcode', 'protein', 'value']
#     df_['type'] = n
#
#     plot_df = concat([plot_df, df_])
#
# plot_df['cnv'] = [cnv_dict[b[:15]][p] if b[:15] in cnv_dict and p in cnv_dict[b[:15]] else np.nan for b, p in plot_df[['barcode', 'protein']].values]
# plot_df = plot_df.dropna()
# print plot_df
#
# sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
# g = sns.FacetGrid(plot_df, row='type', legend_out=True, aspect=1, size=1.5, sharex=True, sharey=True)
# g.map(sns.violinplot, 'value', 'cnv', palette=palette_cnv_number, sym='', linewidth=1., split=False, orient='h')
# g.map(plt.axvline, x=0, ls='--', lw=.3, c='gray')
# g.add_legend()
# g.despine(trim=True)
# g.set_axis_labels('', 'Copy number\nvariation')
# g.set_titles(row_template='{row_name}')
# g.fig.subplots_adjust(wspace=.1, hspace=.4)
# plt.savefig('%s/reports/pancan_cnv_boxplot.pdf' % wd, bbox_inches='tight')
# plt.close('all')
# print '[INFO] Done'

# --
plot_df = DataFrame([])
for n, df in [
    ('residuals', p_residuals),
    ('predicted', p_predicted),
    ('proteomics', proteomics.ix[p_residuals.index]),
    ('transcriptomics', transcriptomics.ix[p_residuals.index].dropna())
]:
    df_ = DataFrame({i: gkn(df.ix[i].dropna()).to_dict() for i in df.index}).T
    # df_ = zscore(df.T).T
    df_ = df_.unstack().reset_index().dropna()
    df_.columns = ['barcode', 'protein', 'value']
    df_['type'] = n

    plot_df = concat([plot_df, df_])

plot_df['cnv'] = [cnv_dict[b[:15]][p] if b[:15] in cnv_dict and p in cnv_dict[b[:15]] else np.nan for b, p in plot_df[['barcode', 'protein']].values]
plot_df = plot_df.dropna()
print plot_df

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.FacetGrid(plot_df, row='type', legend_out=True, aspect=1, size=1.5, sharex=True, sharey=True)
g.map(sns.boxplot, 'value', 'cnv', palette=palette_cnv_number, sym='', linewidth=1., orient='h', notch=True)
g.map(plt.axvline, x=0, ls='--', lw=.3, c='gray')
g.add_legend()
g.despine(trim=True)
g.set_axis_labels('', 'Copy number\nvariation')
g.set_titles(row_template='{row_name}')
g.fig.subplots_adjust(wspace=.1, hspace=.4)
plt.savefig('%s/reports/pancan_cnv_boxplot_residuals.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# --
plot_df = p_residuals.corr()

cmap = sns.diverging_palette(220, 10, n=9, as_cmap=True)
sns.set(style='white', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=0.75)
g = sns.clustermap(
    plot_df, metric='correlation', col_colors=[palette[annot.ix[i, 'type']] for i in plot_df],
    row_colors=[palette[annot.ix[i, 'type']] for i in plot_df.index], cmap=cmap, figsize=(5, 5),
    xticklabels=False, yticklabels=False
)
plt.suptitle('Pancancer proteomics clustermap')
plt.savefig('%s/reports/pancan_clustermap_residuals.png' % wd, bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Plot done'


# --
plot_df = p_residuals.unstack().reset_index().dropna()
plot_df.columns = ['barcode', 'protein', 'residual']
plot_df['type'] = ['Oncogene' if p in cancer_genes['Activating'] else ('Tumor suppressor' if p in cancer_genes['Loss of function'] else 'Neutral') for p in plot_df['protein']]

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.violinplot('residual', 'type', data=plot_df, order=['Neutral', 'Oncogene', 'Tumor suppressor'], linewidth=1., orient='h')
plt.axvline(x=0, ls='--', lw=.3, c='gray')
sns.despine(trim=True)
plt.gcf().set_size_inches(5, 2)
plt.savefig('%s/reports/pancan_residuals_cancer_genes_boxplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# --
protein = 'JAK1'
plot_df = DataFrame({'residual': p_residuals.ix[protein], 'protein': proteomics.ix[protein]})
plot_df['cnv'] = [cnv.ix[protein, i[:15]] if i[:15] in set(cnv) else np.nan for i in plot_df.index]

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.jointplot('cnv', 'residual', data=plot_df, color=default_color, joint_kws={'s': 40, 'edgecolor': 'w', 'linewidth': .5}, kind='scatter')
plt.gcf().set_size_inches(3, 3)
plt.savefig('%s/reports/pancan_residual_plot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
