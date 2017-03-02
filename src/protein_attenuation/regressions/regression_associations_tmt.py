#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv
from statsmodels.stats.multitest import multipletests
from matplotlib_venn import venn2, venn2_circles
from protein_attenuation import palette
from protein_attenuation.utils import read_uniprot_genename

# -- Import phospho
phospho = read_csv('./data/residuals_phospho_protein.csv', index_col=0)


# -- Import gene map
gmap = read_uniprot_genename()

# -- Improt regression results
ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all_tmt.csv')
ppairs_cnv_set = {(px, py) for px, py, fdr in ppairs_cnv[['px', 'py', 'fdr']].values if fdr < .05}

ppairs_phospho = read_csv('./tables/ppairs_phospho_regulation_all_tmt.csv')
ppairs_phospho = ppairs_phospho[[(px.split('_')[0], py) in ppairs_cnv_set for px, py in ppairs_phospho[['px', 'py']].values]]
ppairs_phospho['fdr'] = multipletests(ppairs_phospho['f_pval'], method='fdr_bh')[1]
ppairs_phospho_set = {(px.split('_')[0], py) for px, py, fdr in ppairs_phospho[['px', 'py', 'fdr']].values if fdr < .05}


# --
signif_assoc = ppairs_phospho[[((px.split('_')[0], py) in ppairs_phospho_set) and (fdr < .05) for px, py, fdr in ppairs_phospho[['px', 'py', 'fdr']].values]]
signif_assoc.sort('fdr').to_csv('./tables/significant_phospho_associations.csv', index=False)


# -- Interfaces
interfaces = read_csv('./files/interfaces.tab', sep='\t')
interfaces = interfaces[[s in gmap and t in gmap for s, t in interfaces[['PROTEIN', 'PARTNER']].values]]

interfaces['protein_gene'] = [gmap[p][0] for p in interfaces['PROTEIN']]
interfaces['partner_gene'] = [gmap[p][0] for p in interfaces['PARTNER']]

interfaces['site'] = interfaces['protein_gene'] + '_' + interfaces['AA'] + interfaces['POS'].astype(str)

# px, py = 'MRE11A_S689', 'RAD50'
for px, py in signif_assoc[['px', 'py']].values:
    print interfaces[(interfaces['site'] == px)]
