import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv


# -- Import proteomics
brca_proteomics = read_csv('./data/brca_cell_lines_proteomics_preprocessed.csv', index_col=0)
print brca_proteomics

# -- Signif Px -> Py associations
ppairs = read_csv('./tables/ppairs_cnv_regulation.csv')
ppairs = ppairs.set_index(['px', 'py'])
print ppairs

# --
proteins, samples = set(brca_proteomics.index), set(brca_proteomics)
print 'proteins', 'samples', len(proteins), len(samples)


# --
def pp_correlation(px, py):
    x, y = brca_proteomics.ix[[px, py], samples].dropna(axis=1).values
    cor, pval = pearsonr(x, y)
    return {'px': px, 'py': py, 'cor': cor, 'pval': pval, 'len': len(y)}

ppairs_brca = DataFrame([pp_correlation(px, py) for px, py in ppairs.index if px in proteins and py in proteins])
ppairs_brca['fdr'] = multipletests(ppairs_brca['pval'], method='fdr_bh')[1]
ppairs_brca['ppair_cor'] = [ppairs.ix[(px, py), 'cor'] for px, py in ppairs_brca[['px', 'py']].values]
print ppairs_brca.sort('fdr')


# --
selected_lines = brca_proteomics.ix[set(ppairs_brca['px']).union(ppairs_brca['py']), samples].dropna(axis=1).median().sort_values(ascending=False)
print selected_lines