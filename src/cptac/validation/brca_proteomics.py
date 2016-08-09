import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, default_color, palette_cnv_number, palette
from cptac.utils import gkn
from matplotlib_venn import venn3, venn3_circles
from sklearn.linear_model import LinearRegression
from pandas import DataFrame, Series, read_csv, pivot_table
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Uniprot
uniprot = read_uniprot_genename()

# -- BRCA proteomics
proteomics = read_csv('%s/data/brca_cell_lines_proteomics.csv' % wd)

# Drop tumour samples
proteomics = proteomics[[c for c in proteomics if 'Tumor' not in c]]

# Parse uniprot IDs
proteomics['Uniprot Entry'] = [i.split('|')[1] for i in proteomics['Uniprot Entry']]
proteomics = proteomics[[i in uniprot for i in proteomics['Uniprot Entry']]]
proteomics['Uniprot Entry'] = [uniprot[i][0] for i in proteomics['Uniprot Entry']]
proteomics = proteomics.groupby('Uniprot Entry').mean()

# Log2
proteomics = np.log2(proteomics)

# Average replicates
proteomics.columns = [i.split(' ')[0] for i in proteomics]
proteomics = DataFrame({c: proteomics[c].mean(1) for c in set(proteomics)})

# Proteins measured in at least 50% of the samples
proteomics = proteomics[proteomics.count(1) > proteomics.shape[1] * .5]

# Normalise
proteomics = DataFrame({p: gkn(proteomics.ix[p].dropna()) for p in proteomics.index}).T

# Export
proteomics.to_csv('%s/data/brca_cell_lines_proteomics_preprocessed.csv' % wd)
print 'proteomics', proteomics.shape


# -- GDSC genomics
genomics = read_csv('/Users/emanuel/Projects/data/sanger/genomic/All_variants_cell-lines_30112015_indelMod_ANNOVfun_mutCons.txt', sep='\t')
genomics = genomics[[i in set(proteomics) for i in genomics['SAMPLE_NAME']]]
genomics = genomics[genomics['ZYGOSITY'] == 'homozygous']
genomics['value'] = 1
genomics_m = pivot_table(genomics, index='GENE_NAME', columns='SAMPLE_NAME', values='value', aggfunc='first', fill_value=0)
print 'genomics', genomics_m.shape


# -- Copy number variation
cnv = read_csv('%s/data/CosmicCLP_CompleteCNA.tsv' % wd, sep='\t')
cnv = cnv[[i in set(proteomics) for i in cnv['SAMPLE_NAME']]]
cnv['gene_name'] = [i.split('_')[0] for i in cnv['gene_name']]
cnv['consequence'] = [1 if i == 'gain' else 0 for i in cnv['MUT_TYPE']]
cnv_m = pivot_table(cnv, index='gene_name', columns='SAMPLE_NAME', values='consequence', fill_value=0, aggfunc='first')
print 'genomics', cnv_m.shape


# -- Transcriptomics
samplesheet = read_csv('/Users/emanuel/Projects/data/sanger/samplesheet_paper.tab', index_col=1, sep='\t')

transcriptomics = read_csv('/Users/emanuel/Projects/data/sanger/expression/01_RMAproc_basal_exp.csv', index_col=0)
transcriptomics = transcriptomics.loc[:, [int(i.split('.')[0]) in samplesheet.index for i in transcriptomics]]
transcriptomics.columns = [samplesheet.ix[int(i.split('.')[0]), 'SampleName'] for i in transcriptomics]

transcriptomics = transcriptomics[[i for i in transcriptomics if i in proteomics.columns]]

transcriptomics = DataFrame({g: gkn(transcriptomics.ix[g].dropna()) for g in transcriptomics.index}).T
print 'transcriptomics', transcriptomics.shape


# # Venn
# plot_df = {
#     'CNV': set(cnv_m),
#     'Transcriptomics': set(transcriptomics),
#     'Proteomics': set(proteomics)
# }
#
# venn3(plot_df.values(), set_labels=plot_df.keys(), set_colors=[palette[k] for k in plot_df])
# venn3_circles(plot_df.values(), linestyle='solid', color='white')
# plt.title('GDSC cell lines samples overlap')
# plt.savefig('%s/reports/venn_cell_lines_samples.pdf' % wd, bbox_inches='tight')
# plt.close('all')
# print '[INFO] Done'


# -- Overlap
samples = set(proteomics).intersection(transcriptomics)
proteins = set(proteomics.index).intersection(transcriptomics.index)
print 'samples', 'proteins', len(samples), len(proteins)

proteomics, transcriptomics = proteomics.ix[proteins, samples], transcriptomics.ix[proteins, samples]
print 'proteomics', 'transcriptomics', proteomics.shape, transcriptomics.shape


# -- Residual
def protein_residual(g):
    y = proteomics.ix[g, samples].dropna()
    x = transcriptomics.ix[[g], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[g] - lm.intercept_

    return y_
residuals = DataFrame({g: protein_residual(g) for g in proteins}).T
print 'residuals', residuals.shape


# -- Protein-pairs
ppairs = read_csv('%s/tables/ppairs_cnv_regulation.csv' % wd)
ppairs = ppairs[ppairs['cor'] > 0]
ppairs = ppairs[[px in proteins and py in proteins for px, py in ppairs[['px', 'py']].values]]
print ppairs


# --
plot_df = DataFrame([{'px': px, 'py': py, 'sample': s, 'x': residuals.ix[px, s], 'y': residuals.ix[py, s]} for px, py in ppairs[['px', 'py']].values for s in proteomics]).dropna()
plot_df['mutation_x'] = [cnv_m.ix[p, s] if p in cnv_m.index else 0 for p, s in plot_df[['px', 'sample']].values]
print plot_df.sort(['mutation_x', 'y'])

sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
sns.boxplot('mutation_x', 'y', data=plot_df, sym='', linewidth=.3, palette=palette_cnv_number)
sns.stripplot('mutation_x', 'y', data=plot_df, split=True, jitter=True, edgecolor='white', linewidth=.3, palette=palette_cnv_number)
plt.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
sns.despine(trim=True)
plt.xlabel('Copy number')
plt.ylabel('Residuals')
plt.gcf().set_size_inches(1, 4)
plt.savefig('%s/reports/ppairs_cnv_regulation_validation_brca_boxplots.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
