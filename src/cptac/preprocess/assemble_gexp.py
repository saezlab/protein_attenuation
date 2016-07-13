import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from pandas import DataFrame, Series, read_csv

prot = {i[:16] for i in read_csv('%s/tables/pancan_preprocessed_normalised.csv' % wd, index_col=0)}
print prot

gexp = read_csv('/Users/emanuel/Projects/data/cptac/GSE62944_merged_expression_voom.tsv', sep='\t', index_col=0)
gexp = gexp.loc[:, [i[:16] in prot for i in gexp]]
gexp.to_csv('%s/data/tcga_rnaseq.tsv' % wd, sep='\t')
print gexp



