from cptac import wd
from pandas import DataFrame, Series, read_csv

# Overlapping samples with the proteome
samples = set(Series.from_csv('%s/data/samplesheet.csv' % wd).index)
print len(samples)


# -- Voom transformed
# Import TCGA pancancer rna-seq data
gexp = read_csv('/Users/emanuel/Projects/data/cptac/GSE62944_merged_expression_voom.tsv', sep='\t', index_col=0)
print gexp.shape

# Overlap
gexp = gexp.loc[:, [i[:15] in samples for i in gexp]]

# Discard metastatic
gexp = gexp[[i for i in gexp if i[15] == 'A']]
gexp.columns = [c[:15] for c in gexp]

# Average replicates
gexp = DataFrame({i: gexp.loc[:, [i]].mean(1) for i in set(gexp)})

# Export
gexp.to_csv('%s/data/tcga_rnaseq.tsv' % wd, sep='\t')
print gexp


# -- FPKM
fpkm = read_csv('%s/files/GSM1536837_01_27_15_TCGA_20.Illumina.tumor_Rsubread_FPKM.txt' % wd, sep='\t', index_col=0)
print 'fpkm', fpkm.shape

# Overlap
fpkm = fpkm.loc[:, [i[:15] in samples for i in fpkm]]

# Discard metastatic
fpkm = fpkm[[i for i in fpkm if i[15] == 'A']]
fpkm.columns = [c[:15] for c in fpkm]

# Average replicates
fpkm = DataFrame({i: fpkm.loc[:, [i]].mean(1) for i in set(fpkm)})

# Export
fpkm.to_csv('%s/data/tcga_rnaseq_fpkm.csv' % wd)
print fpkm
