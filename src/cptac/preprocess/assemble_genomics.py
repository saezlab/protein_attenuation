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

# Overlapping samples with the proteome
samples = set(Series.from_csv('%s/data/samplesheet.csv' % wd).index)
print len(samples)


# Source data folder
dwd = '/Users/emanuel/Projects/data/cptac/tcga_mutations_firehose_06_06_2016_maf_files/'

# Read existing maf files
mafs = [i for i in os.listdir(dwd) if i.endswith('.maf.txt')]
print len(mafs)

# Merge maf files
genomics = concat([read_csv('%s/%s' % (dwd, f), sep='\t') for f in mafs]).reset_index()
print len(set(genomics['Tumor_Sample_Barcode']))

# Sub-set by proteomics samples
genomics = genomics[[i[:15] in samples for i in genomics['Tumor_Sample_Barcode']]]
print len({i[:15] for i in genomics['Tumor_Sample_Barcode']})

# Export cnv
genomics.to_csv('%s/data/tcga_genomics.tsv' % wd, sep='\t')
print '[INFO] Exported'
