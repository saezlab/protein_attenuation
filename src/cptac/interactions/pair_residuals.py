import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from cptac import wd, genomic_mod, palette
from sklearn.linear_model import LinearRegression
from pandas import DataFrame, Series, read_csv, pivot_table, concat


# -- Import top correlated proteins
pip = read_csv('%s/tables/top_correlated_protein_pairs.csv' % wd)
print pip

# -- Create network
network_i = igraph.Graph(directed=False)

# Add nodes
network_i.add_vertices(list({i for p in pip[['p1', 'p2']].values for i in p}))

# Add edges
network_i.add_edges([(p1, p2) for p1, p2 in pip[['p1', 'p2']].values])

# Remove duplicated edges and self-loops
network_i = network_i.simplify(True, True, 'first')
print Series(dict(zip(*(network_i.vs['name'], network_i.degree())))).sort_values()

# -- Import genomics
genomics = read_csv('%s/tables/tcga_genomics.csv' % wd)
genomics = genomics[[i in genomic_mod for i in genomics['Variant_Classification']]]
genomics['patient_id'] = ['-'.join(i.split('-')[:4])[:-1] for i in genomics['Tumor_Sample_Barcode']]
genomics['value'] = 1
print genomics

# -- Import pancancer proteomics
pancan = read_csv('%s/tables/pancan_preprocessed_normalised.csv' % wd, index_col=0)
pancan = pancan[pancan.count(1) > (pancan.shape[1] * .5)]
pancan_samplesheet = read_csv('%s/tables/samplesheet.csv' % wd, index_col=0)
print pancan


# -- Estimate residuals
def residuals(p1, p2):
    name = '%s-%s' % (p1, p2)
    print name

    # Get protein measurements
    x = pancan.ix[[p1, p2]].T.dropna()[p1]
    y = pancan.ix[[p1, p2]].T.dropna()[p2]

    # Add tissue covariates
    x = concat([sm.add_constant(x), pancan_samplesheet.ix[y.index, 'type'].str.get_dummies()], axis=1)

    # Fit linear regression
    lm = sm.OLS(y, x).fit()

    return lm.resid

p_residuals = DataFrame({'%s_%s' % (p1, p2): residuals(p1, p2) for p1, p2 in pip[['p1', 'p2']].values}).T.round(5)
p_residuals.to_csv('%s/tables/protein_pairs_residuals.csv' % wd)
print p_residuals


#
tot_mut = pivot_table(genomics, index='Hugo_Symbol', columns='patient_id', values='value', aggfunc=np.sum, fill_value=0)

var_df = DataFrame({'residuals': np.log10(p_residuals.abs().sum())})
var_df['mutations'] = [np.log10(tot_mut.sum().ix['-'.join(i.split('-')[:4])[:-1]]) if '-'.join(i.split('-')[:4])[:-1] in tot_mut.columns else np.nan for i in var_df.index]
var_df['type'] = [pancan_samplesheet.ix[i, 'type'] for i in var_df.index]
var_df = var_df.dropna()
print var_df.corr()

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=.75)
for g in set(var_df['type']):
    sns.regplot(var_df.ix[var_df['type'] == g, 'residuals'], var_df.ix[var_df['type'] == g, 'mutations'], label=g, color=palette[g], fit_reg=False)
sns.despine()

plt.legend()
plt.savefig('%s/reports/protein_pair_residuals_genomics_scatter.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] done'


#
def rm_batch(x, y):
    ys = y.dropna()
    xs = x.ix[ys.index]

    lm = LinearRegression().fit(xs.values, ys.values)

    return ys - xs.dot(lm.coef_) - lm.intercept_

design = var_df['type'].str.get_dummies()
var_df_n = DataFrame({p: rm_batch(design, var_df.ix[design.index, p]) for p in ['residuals', 'mutations']})
var_df_n['type'] = [pancan_samplesheet.ix[i, 'type'] for i in var_df_n.index]
print var_df_n.corr()


sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=.75)
for g in set(var_df_n['type']):
    sns.regplot(var_df_n.ix[var_df_n['type'] == g, 'residuals'], var_df_n.ix[var_df_n['type'] == g, 'mutations'], label=g, color=palette[g], fit_reg=False)
sns.despine()

plt.legend()
plt.savefig('%s/reports/protein_pair_residuals_genomics_scatter_no_tissue.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] done'
