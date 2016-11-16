#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import matplotlib.pyplot as plt
from pandas import read_csv
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict

# -- Import CORUM
uniprot = read_uniprot_genename()

corum_n = get_complexes_name()

# dict
corum = {k: {uniprot[p][0] for p in v if p in uniprot} for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if 1 < len(corum[k])}
print 'corum', len(corum)

# -- Import RPPA list
rppa = read_csv('./tables/rppa_ab_gene_info.csv')
rppa['corum_id'] = [';'.join([str(c) for c in corum if i in corum[c]]) for i in rppa['Gene']]
rppa['corum_name'] = [';'.join([corum_n[c] for c in corum if i in corum[c]]) for i in rppa['Gene']]
rppa.to_csv('./tables/rppa_ab_gene_info_annotated.csv', index=False)
print rppa
