import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from cptac import wd, palette
from pandas import DataFrame, Series, read_csv, pivot_table

# Overlapping samples with the proteome
samples = set(Series.from_csv('%s/data/samplesheet.csv' % wd).index)
print len(samples)

# Import whole CNV data-set
cnv = read_csv('/Users/emanuel/Projects/data/cptac/cna_thresholded.tsv', sep='\t')
print cnv

gexp = read_csv('/Users/emanuel/Projects/data/cptac/GSE62944_merged_expression_voom.tsv', sep='\t', index_col=0)
print gexp.shape


# Venn diagram
plot_df = {
    'CNV': {i[:15] for i in cnv['barcode']},
    'Transcriptomics': {i[:15] for i in gexp},
    'Proteomics': samples
}

venn3(plot_df.values(), set_labels=plot_df.keys(), set_colors=[palette[k] for k in plot_df])
venn3_circles(plot_df.values(), linestyle='solid', color='white')
plt.title('Overlap between proteomics and copy number variation')
plt.savefig('%s/reports/venn_cnv_proteomics.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'

# Sub-set by proteomics samples
cnv = cnv[[i[:15] in samples for i in cnv['barcode']]]
print len({i[:15] for i in cnv['barcode']})

# Build matrix - duplicated entries on same sample are discarded if gistic values differ
cnv = cnv.groupby(['barcode', 'hgnc'])['gistic'].agg(lambda x: np.nan if len(set(x)) > 1 else list(x)[0])
cnv = cnv.reset_index()
cnv = cnv.dropna()
cnv = pivot_table(cnv, index='hgnc', columns='barcode', values='gistic', fill_value=0)
print cnv

# Parse sample ID
cnv.columns = [i[:15] for i in cnv]

# Export cnv
cnv.to_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t')
print '[INFO] Exported'
