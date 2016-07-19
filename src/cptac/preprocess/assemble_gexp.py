from cptac import wd
from pandas import DataFrame, Series, read_csv


# Overlapping samples with the proteome
samples = set(Series.from_csv('%s/data/samplesheet.csv' % wd).index)
print len(samples)

# Import TCGA pancancer rna-seq data
gexp = read_csv('/Users/emanuel/Projects/data/cptac/GSE62944_merged_expression_voom.tsv', sep='\t', index_col=0)
print gexp.shape

# Overlap
gexp = gexp.loc[:, [i[:15] in samples for i in gexp]]
gexp.columns = [c[:15] for c in gexp]

# Average replicates
gexp = DataFrame({i: gexp.loc[:, [i]].mean(1) for i in set(gexp)})

# Export
gexp.to_csv('%s/data/tcga_rnaseq.tsv' % wd, sep='\t')
print gexp
