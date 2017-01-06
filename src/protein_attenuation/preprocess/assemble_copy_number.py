#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
from pandas import Series, read_csv, pivot_table

# Overlapping samples with the proteome
samples = set(Series.from_csv('./data/samplesheet.csv').index)

# Import whole CNV data-set
cnv = read_csv('./data/tcga_pancan_cna_gistic_thresholded.tsv', sep='\t')

# Consider only primary tumour samples
cnv = cnv[[i[13:16] == '01A' for i in cnv['barcode']]]

# Sub-set by proteomics samples
cnv = cnv[[i[:12] in samples for i in cnv['barcode']]]

# Build matrix - duplicated entries on same sample are discarded if gistic values differ
cnv = cnv.groupby(['barcode', 'hgnc'])['gistic'].agg(lambda x: np.nan if len(set(x)) > 1 else list(x)[0])
cnv = cnv.reset_index()
cnv = cnv.dropna()
cnv = pivot_table(cnv, index='hgnc', columns='barcode', values='gistic', fill_value=0)

# Parse sample ID
cnv.columns = [i[:12] for i in cnv]

# Export cnv
cnv.to_csv('./data/tcga_cnv.csv')
print '[INFO] Copy-number exported: ' + './data/tcga_cnv.csv'
