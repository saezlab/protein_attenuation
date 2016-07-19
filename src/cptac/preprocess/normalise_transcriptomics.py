from cptac import wd
from cptac.utils import gkn
from pandas import DataFrame, read_csv


# -- Normalise transcriptomics
for f in ['tcga_rnaseq', 'tcga_rnaseq_corrected']:
    transcriptomics = read_csv('%s/data/%s.tsv' % (wd, f), sep='\t', index_col=0)
    transcriptomics = DataFrame({i: gkn(transcriptomics.ix[i].dropna()).to_dict() for i in transcriptomics.index}).T
    transcriptomics.to_csv('%s/data/%s_normalised.csv' % (wd, f))
    print '[INFO] Done: %s' % f
