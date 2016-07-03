import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from pandas import DataFrame, Series, read_csv


# -- BRCA
brca = read_csv('%s/data/brca_proteomics.csv' % wd)

# Remove ambigously mapping peptides
brca = brca[[len(i.split('|')) == 1 for i in brca['accession_numbers']]]

# Remove non-gene mapping proteins and average duplicated
brca = brca.dropna(subset=['geneName']).groupby('geneName').mean()

# Match sample ids
brca_samplesheet = read_csv('%s/data/brca_proteomics_samplesheet.csv' % wd, index_col=0)

brca = brca.loc[:, [i.split('.')[0] in brca_samplesheet.index for i in brca]]
brca.columns = [brca_samplesheet.ix[i.split('.')[0], 'Biospecimen Barcode Side'] for i in brca]

# Remove proteins measured in less than 25% of the samples
brca = brca[(brca.count(1) / brca.shape[1]) > .25]

# Export
brca.to_csv('%s/tables/brca_proteomics_processed.csv' % wd)
print '[INFO] BRCA proteomics preprocessed exported'


# -- HGSC
hgsc = read_csv('%s/data/hgsc_proteomics.csv' % wd).drop('refseq_peptide', axis=1)

# Remove non-gene mapping proteins and average duplicated
hgsc = hgsc.dropna(subset=['hgnc_symbol']).groupby('hgnc_symbol').mean()

# Match sample ids
hgsc_samplesheet = read_csv('%s/data/hgsc_proteomics_samplesheet.csv' % wd)
hgsc_samplesheet['id'] = ['%s-%s' % (i, c) for i, c in hgsc_samplesheet[['PCC', 'bcr_patient_barcode']].values]
hgsc_samplesheet = hgsc_samplesheet.set_index('id')

hgsc.columns = [hgsc_samplesheet.ix[i, 'TCGA barcode (shipped portion)'] for i in hgsc]

# Remove proteins measured in less than 25% of the samples
hgsc = hgsc[(hgsc.count(1) / hgsc.shape[1]) > .25]

# Export
hgsc.to_csv('%s/tables/hgsc_proteomics_processed.csv' % wd)
print '[INFO] HGSC proteomics preprocessed exported'


# -- COREAD
coread = read_csv('%s/data/coread_proteomics.csv' % wd, index_col=0).replace(0, np.nan)

# Remove proteins measured in less than 25% of the samples
coread = coread[(coread.count(1) / coread.shape[1]) > .25]

# Export
coread.to_csv('%s/tables/coread_proteomics_processed.csv' % wd)
print '[INFO] COREAD proteomics preprocessed exported'
