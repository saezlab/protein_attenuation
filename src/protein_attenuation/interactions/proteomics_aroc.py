#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import itertools as it
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from protein_attenuation import wd, palette, palette_dbs
from gdsc import wd as gdsc_wd
from matplotlib.gridspec import GridSpec
from pandas import read_csv, DataFrame, concat, Series
from sklearn.metrics.classification import matthews_corrcoef
from sklearn.metrics.ranking import roc_curve, auc, precision_recall_curve, average_precision_score
from pymist.utils.stringdb import get_stringdb, get_stringdb_actions
from pymist.utils.biogriddb import get_biogriddb, get_biogriddb_action
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
uniprot = read_uniprot_genename()

# CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print 'corum', len(corum)

# String
string_thres = {'low': 150, 'medium': 400, 'high': 700, 'highest': 900}

string = {thres: {(uniprot[s][0], uniprot[t][0]) for p1, p2 in get_stringdb(string_thres[thres]) for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot} for thres in string_thres}
print 'string', len(string)

string_action = {thres: get_stringdb_actions(string_thres[thres]).groupby('mode')['interaction'].agg(lambda x: set(x)).to_dict() for thres in string_thres}
string_action = {thres: {mode: {(uniprot[s][0], uniprot[t][0]) for p1, p2 in string_action[thres][mode] for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot} for mode in string_action[thres]} for thres in string_thres}
print 'string_action', len(string_action)

# Biogrid
biogrid = get_biogriddb()
biogrid = {(s, t) for p1, p2 in biogrid for s, t in [(p1, p2), (p2, p1)]}
print 'biogrid', len(biogrid)

biogrid_action = get_biogriddb_action().groupby('Experimental System')['interaction'].agg(lambda x: set(x)).to_dict()
print 'biogrid_action', len(biogrid_action)

# Omnipath
omnipath = read_csv('%s/files/omnipath_interactions.txt' % wd, sep='\t')
omnipath = {(uniprot[s][0], uniprot[t][0]) for s, t in omnipath[['source', 'target']].values if s in uniprot and t in uniprot}
print 'omnipath', len(omnipath)

omnipath_action = read_csv('%s/files/omnipath_interactions.txt' % wd, sep='\t')
omnipath_sources = {s for i in omnipath_action['sources'] for s in i.split(';')}
omnipath_action = {s: {(uniprot[s][0], uniprot[t][0]) for s, t in omnipath_action[[s in i for i in omnipath_action['sources']]][['source', 'target']].values if s in uniprot and t in uniprot} for s in omnipath_sources}
print 'omnipath_action', len(omnipath_action)


# -- Proteomics
brca = read_csv('%s/data/brca_proteomics_processed.csv' % wd, index_col=0)
brca = brca[brca.count(1) > (brca.shape[1] * .5)]

hgsc = read_csv('%s/data/hgsc_proteomics_processed.csv' % wd, index_col=0)
hgsc = hgsc[hgsc.count(1) > (hgsc.shape[1] * .5)]

coread = read_csv('%s/data/coread_proteomics_processed.csv' % wd, index_col=0)
coread = coread[coread.count(1) > (coread.shape[1] * .5)]

ov_prot = set(brca.index).intersection(hgsc.index).intersection(coread.index)
brca, hgsc, coread = brca.ix[ov_prot], hgsc.ix[ov_prot], coread.ix[ov_prot]

proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0).ix[ov_prot]
proteomics = proteomics[proteomics.count(1) > (proteomics.shape[1] * .5)]
print 'proteomics', proteomics.shape

transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0).ix[ov_prot].dropna()
print 'transcriptomics', transcriptomics.shape


# -- Protein-protein correlation
cor_res, cor_dfs = {}, {}
# name, d_df = 'Proteomics', proteomics
for name, d_df in [('BRCA', brca), ('HGSC', hgsc), ('COREAD', coread), ('Proteomics', proteomics), ('Transcriptomics', transcriptomics)]:
    print name

    # -- Generate p-p correlations
    df = d_df.ix[{p for x, y in corum for p in [x, y]}.intersection(proteomics.index)].T.corr(method='pearson')
    df.values[np.tril_indices(df.shape[0], 0)] = np.nan
    df.index.name = None
    df = df.unstack().reset_index().dropna()
    df.columns = ['p1', 'p2', 'cor']

    # -- Annotate p-p interactions
    df['CORUM'] = [1 if ((p1, p2) in corum) else 0 for p1, p2 in df[['p1', 'p2']].values]
    df['STRING'] = [1 if ((p1, p2) in string['highest']) else 0 for p1, p2 in df[['p1', 'p2']].values]
    df['BioGRID'] = [1 if ((p1, p2) in biogrid) else 0 for p1, p2 in df[['p1', 'p2']].values]
    df['OmniPath'] = [1 if ((p1, p2) in omnipath) else 0 for p1, p2 in df[['p1', 'p2']].values]

    for thres in string_thres:
        df['STRING_%s' % thres] = [1 if ((p1, p2) in string[thres]) else 0 for p1, p2 in df[['p1', 'p2']].values]
        print 'STRING_%s' % thres

    for thres in string_action:
        for action in string_action[thres]:
            df['STRING_%s_%s' % (thres, action)] = [1 if ((p1, p2) in string_action[thres][action]) else 0 for p1, p2 in df[['p1', 'p2']].values]
            print 'STRING_%s_%s' % (thres, action)

    for action in biogrid_action:
        df['BioGRID_%s' % action] = [1 if ((p1, p2) in biogrid_action[action]) else 0 for p1, p2 in df[['p1', 'p2']].values]
        print 'BioGRID_%s' % action

    for action in omnipath_action:
        df['OmniPath_%s' % action] = [1 if ((p1, p2) in omnipath_action[action]) else 0 for p1, p2 in df[['p1', 'p2']].values]
        print 'OmniPath_%s' % action

    # -- Estimate ROC curves
    cor_res[name] = {}
    for db in ['CORUM', 'STRING', 'BioGRID', 'OmniPath']:
        curve_fpr, curve_tpr, _ = roc_curve(df[db], df['cor'])
        cor_res[name][db] = (curve_fpr, curve_tpr)

    for thres in string_thres:
        curve_fpr, curve_tpr, _ = roc_curve(df['STRING_%s' % thres], df['cor'])
        cor_res[name]['STRING_%s' % thres] = (curve_fpr, curve_tpr)

    for thres in string_action:
        for action in string_action[thres]:
            curve_fpr, curve_tpr, _ = roc_curve(df['STRING_%s_%s' % (thres, action)], df['cor'])
            cor_res[name]['STRING_%s_%s' % (thres, action)] = (curve_fpr, curve_tpr)

    for action in biogrid_action:
        curve_fpr, curve_tpr, _ = roc_curve(df['BioGRID_%s' % action], df['cor'])
        cor_res[name]['BioGRID_%s' % action] = (curve_fpr, curve_tpr)

    for action in omnipath_action:
        curve_fpr, curve_tpr, _ = roc_curve(df['OmniPath_%s' % action], df['cor'])
        cor_res[name]['OmniPath_%s' % action] = (curve_fpr, curve_tpr)

    # -- Store data-frame
    cor_dfs[name] = df

    # -- Histograms
    sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'}, font_scale=.75)
    fig, gs = plt.figure(figsize=(20, 3)), GridSpec(1, 5, hspace=.3, wspace=.3)

    # All data-bases
    ax = plt.subplot(gs[0])
    g = sns.distplot(df['cor'], color=palette_dbs['All'], hist=False, kde_kws={'shade': True}, label='All', ax=ax)
    for db in ['CORUM', 'STRING', 'BioGRID', 'OmniPath']:
        sns.distplot(df.ix[df[db] == 1, 'cor'], color=palette_dbs[db], hist=False, kde_kws={'shade': True}, label=db, ax=ax)

    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_title('%s\ndata-bases' % name)
    ax.set_xlabel('Pearson\'s r')
    ax.set_ylabel('Density')
    g.set_xlim(-1, 1)
    sns.despine(trim=True, ax=ax)

    # STRING thresholds
    ax = plt.subplot(gs[1])
    order = ['low', 'medium', 'high', 'highest']
    pal = dict(zip(*(order, sns.light_palette(palette_dbs['STRING'], n_colors=len(order) + 1).as_hex()[1:])))
    for thres in order:
        g = sns.distplot(df.ix[df['STRING_%s' % thres] == 1, 'cor'], color=pal[thres], hist=False, kde_kws={'shade': True}, label='%s (> %d)' % (thres, string_thres[thres]), ax=ax)
        g.set_xlim(-1, 1)

    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_title('%s\nSTRING confidence thresholds' % name)
    ax.set_xlabel('Pearson\'s r')
    ax.set_ylabel('Density')
    sns.despine(trim=True, ax=ax)

    # STRING type
    ax = plt.subplot(gs[2])
    thres = 'highest'
    order = ['reaction', 'binding', 'catalysis', 'activation']
    pal = dict(zip(*(order, sns.light_palette(palette_dbs['STRING'], n_colors=len(order) + 1).as_hex()[1:])))
    for action in order:
        if len(df.ix[df['STRING_%s_%s' % (thres, action)] == 1, 'cor']) > 0:
            g = sns.distplot(df.ix[df['STRING_%s_%s' % (thres, action)] == 1, 'cor'], color=pal[action], hist=False, kde_kws={'shade': True}, label=action, ax=ax)
            g.set_xlim(-1, 1)

    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_title('%s\nSTRING interaction type (threshold: %s)' % (name, thres))
    ax.set_xlabel('Pearson\'s r')
    ax.set_ylabel('Density')
    sns.despine(trim=True, ax=ax)

    # BioGRID
    ax = plt.subplot(gs[3])
    order = ['Co-purification', 'Co-crystal Structure', 'Far Western', 'FRET']
    pal = dict(zip(*(order, sns.light_palette(palette_dbs['BioGRID'], n_colors=len(order) + 1).as_hex()[1:])))
    for action in order:
        if len(df.ix[df['BioGRID_%s' % action] == 1, 'cor']) > 0:
            g = sns.distplot(df.ix[df['BioGRID_%s' % action] == 1, 'cor'], color=pal[action], hist=False, kde_kws={'shade': True}, label=action, ax=ax)
            g.set_xlim(-1, 1)

    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_title('%s\nBioGRID interaction type' % name)
    ax.set_xlabel('Pearson\'s r')
    ax.set_ylabel('Density')
    sns.despine(trim=True, ax=ax)

    # OmniPath
    ax = plt.subplot(gs[4])
    order = ['KEGG', 'PhosphoSite', 'Signor', 'Laudanna_sigflow']
    pal = dict(zip(*(order, sns.light_palette(palette_dbs['OmniPath'], n_colors=len(order) + 1).as_hex()[1:])))
    for action in order:
        if len(df.ix[df['OmniPath_%s' % action] == 1, 'cor']) > 0:
            g = sns.distplot(df.ix[df['OmniPath_%s' % action] == 1, 'cor'], color=pal[action], hist=False, kde_kws={'shade': True}, label=action, ax=ax)
            g.set_xlim(-1, 1)

    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_title('%s\nOmniPath resource' % name)
    ax.set_xlabel('Pearson\'s r')
    ax.set_ylabel('Density')
    sns.despine(trim=True, ax=ax)

    # Export
    plt.savefig('%s/reports/correlation_histogram_%s.pdf' % (wd, name), bbox_inches='tight')
    plt.close('all')

    # -- Export top correlated proteins list
    df[df['cor'] > .4].to_csv('%s/tables/top_correlated_protein_pairs_%s.csv' % (wd, name.lower()), index=False)

    print 'ROC done'

print '[INFO] Done'


# -- Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
fig, gs, pos = plt.figure(figsize=(13, 3)), GridSpec(1, 4, hspace=.3, wspace=.3), 0

for db in ['CORUM', 'STRING', 'BioGRID', 'OmniPath']:
    ax = plt.subplot(gs[pos])

    for name in ['BRCA', 'COREAD', 'HGSC', 'Proteomics', 'Transcriptomics']:
        curve_fpr, curve_tpr = cor_res[name][db]
        plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (name, auc(curve_fpr, curve_tpr)), c=palette[name])

    ax.plot([0, 1], [0, 1], 'k--', lw=.3)
    sns.despine(trim=True)
    ax.legend(loc='lower right')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title(db)

    pos += 1

plt.savefig('%s/reports/proteomics_ppi_aroc.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Plot STRING thresholds
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
fig, gs, pos = plt.figure(figsize=(13, 3)), GridSpec(1, 4, hspace=.3, wspace=.3), 0

for thres in ['highest', 'high', 'medium', 'low']:
    ax = plt.subplot(gs[pos])

    for name in ['BRCA', 'COREAD', 'HGSC', 'Proteomics', 'Transcriptomics']:
        curve_fpr, curve_tpr = cor_res[name]['STRING_%s' % thres]
        plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (name, auc(curve_fpr, curve_tpr)), c=palette[name])

    ax.plot([0, 1], [0, 1], 'k--', lw=.3)
    sns.despine(trim=True)
    ax.legend(loc='lower right')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title('STRING\nthreshold: %s' % thres)

    pos += 1

plt.savefig('%s/reports/proteomics_ppi_aroc_string_threshold.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Plot STRING interactions types by threshold
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs, pos = GridSpec(7, 4, hspace=.45, wspace=.3), 0

for action in ['activation', 'inhibition', 'binding', 'catalysis', 'expression', 'ptmod', 'reaction']:
    for thres in ['highest', 'high', 'medium', 'low']:
        ax = plt.subplot(gs[pos])

        for name in ['BRCA', 'COREAD', 'HGSC', 'Proteomics', 'Transcriptomics']:
            if 'STRING_%s_%s' % (thres, action) in cor_res[name]:
                curve_fpr, curve_tpr = cor_res[name]['STRING_%s_%s' % (thres, action)]
                plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (name, auc(curve_fpr, curve_tpr)), c=palette[name])

        ax.plot([0, 1], [0, 1], 'k--', lw=.3)
        sns.despine(trim=True)
        ax.legend(loc='lower right')
        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.set_title('STRING\nthreshold: %s, type: %s' % (thres, action))

        pos += 1

plt.gcf().set_size_inches(13, 20)
plt.savefig('%s/reports/proteomics_ppi_aroc_string_type.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Plot BioGRID interactions types
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs, pos = GridSpec(5, 5, hspace=.45, wspace=.45), 0

for action in biogrid_action:
    ax = plt.subplot(gs[pos])

    for name in ['BRCA', 'COREAD', 'HGSC', 'Proteomics', 'Transcriptomics']:
        curve_fpr, curve_tpr = cor_res[name]['BioGRID_%s' % action]
        plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (name, auc(curve_fpr, curve_tpr)), c=palette[name])

    ax.plot([0, 1], [0, 1], 'k--', lw=.3)
    sns.despine(trim=True)
    ax.legend(loc='lower right')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title('BioGRID\ntype: %s' % action)

    pos += 1

plt.gcf().set_size_inches(18, 18)
plt.savefig('%s/reports/proteomics_ppi_aroc_biogrid_type.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Plot OmniPath interactions by resources types
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs, pos = GridSpec(8, 5, hspace=.45, wspace=.45), 0

for action in omnipath_action:
    ax = plt.subplot(gs[pos])

    for name in ['BRCA', 'COREAD', 'HGSC', 'Proteomics', 'Transcriptomics']:
        curve_fpr, curve_tpr = cor_res[name]['OmniPath_%s' % action]
        plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (name, auc(curve_fpr, curve_tpr)), c=palette[name])

    ax.plot([0, 1], [0, 1], 'k--', lw=.3)
    sns.despine(trim=True)
    ax.legend(loc='lower right')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title('OmniPath \nresource: %s' % action)

    pos += 1

plt.gcf().set_size_inches(18, 29)
plt.savefig('%s/reports/proteomics_ppi_aroc_omnipath_type.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
