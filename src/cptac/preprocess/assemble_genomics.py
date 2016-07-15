import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, palette_dbs, genomic_mod
from pymist.utils.stringdb import get_stringdb
from pymist.utils.biogriddb import get_biogriddb
from pymist.utils.corumdb import get_complexes_pairs
from pandas import DataFrame, Series, read_csv, pivot_table, concat
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# Source data folder
dwd = '/Users/emanuel/Projects/data/cptac/tcga_mutations_firehose_06_06_2016_maf_files/'

# Read existing maf files
mafs = [i for i in os.listdir(dwd) if i.endswith('.maf.txt')]
print len(mafs)

# Merge maf files
genomics = concat([read_csv('%s/%s' % (dwd, f), sep='\t') for f in mafs]).reset_index()
print len(set(genomics['Tumor_Sample_Barcode']))

# Filter mutation type
genomics = genomics[[i in genomic_mod for i in genomics['type']]]
print genomics

# Overlap with proteomics samples
samples = {i[:15] for i in set(read_csv('%s/data/pancan_preprocessed_normalised.csv' % wd, index_col=0))}

print len(set(genomics['Tumor_Sample_Barcode']).intersection(samples))

# genomics = genomics[[i in samples for i in genomics['Tumor_Sample_Barcode']]]
