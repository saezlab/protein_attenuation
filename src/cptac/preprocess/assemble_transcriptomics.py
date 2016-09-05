#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

from pandas import DataFrame, Series, read_csv


# Overlapping samples with the proteome
samples = set(Series.from_csv('./data/samplesheet.csv').index)
print len(samples)

# -- Voom transformed
# Import TCGA pancancer rna-seq data
gexp = read_csv('/Users/emanuel/Projects/data/cptac/GSE62944_merged_expression_voom.tsv', sep='\t', index_col=0)
print gexp.shape

# Consider only primary tumour samples
gexp = gexp[[i for i in gexp if i[13:16] == '01A']]

# Overlap
gexp = gexp.loc[:, [i[:12] in samples for i in gexp]]
gexp.columns = [c[:12] for c in gexp]

# Average replicates
gexp = DataFrame({i: gexp.loc[:, [i]].mean(1) for i in set(gexp)})

# Export
gexp.to_csv('./data/tcga_rnaseq.csv')
print gexp
