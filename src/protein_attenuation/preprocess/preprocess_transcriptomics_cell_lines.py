#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

from pandas import DataFrame, read_csv

# Samplesheet
samplesheet = read_csv('./data/sanger_samplesheet.csv', index_col=1)

# Gene expression
trans = read_csv('./data/sanger_gene_experssion_rma.txt', sep='\t').dropna().drop('GENE_title', axis=1).set_index('GENE_SYMBOLS')

# Average replicates
trans.columns = ['.'.join(c.split('.')[:-1]) if len(c.split('.')) == 3 else c for c in trans]
trans = DataFrame({i: trans.loc[:, [i]].mean(1) for i in set(trans)})

# Cell lines in samplesheet
trans = trans.loc[:, [int(c.split('.')[1]) in samplesheet.index for c in trans]]

# Map cell lines name
trans.columns = [samplesheet.ix[int(c.split('.')[1]), 'sample'] for c in trans]

# Export
trans.to_csv('./data/sanger_gene_experssion_rma_preprocessed.csv')
print '[INFO] Cell lines microarry gene expression: ', './data/sanger_gene_experssion_rma_preprocessed.csv'
