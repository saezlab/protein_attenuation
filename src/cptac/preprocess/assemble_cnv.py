import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from pandas import DataFrame, Series, read_csv, pivot_table

# Import whole CNV data-set
cnv = read_csv('/Users/emanuel/Projects/data/cptac/cna_thresholded.tsv', sep='\t')
print cnv

# Sub-set by proteomics samples
samples = {i[:15] for i in set(read_csv('%s/tables/pancan_preprocessed_normalised.csv' % wd, index_col=0))}

cnv = cnv[[i[:15] in samples for i in cnv['barcode']]]
print len({i[:15] for i in cnv['barcode']})

# Build matrix - duplicated entries on same sample are discarded if gistic values differ
cnv = cnv.groupby(['barcode', 'hgnc'])['gistic'].agg(lambda x: np.nan if len(set(x)) > 1 else list(x)[0])
cnv = cnv.reset_index()
cnv = cnv.dropna()
cnv = pivot_table(cnv, index='hgnc', columns='barcode', values='gistic', fill_value=0)
print cnv

# Export cnv
cnv.to_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t')
print '[INFO] Exported'
