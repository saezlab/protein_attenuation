#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from protein_attenuation import palette
from matplotlib.gridspec import GridSpec
from pandas import DataFrame, read_csv, concat
from sklearn.metrics.ranking import roc_curve, auc
from protein_attenuation.utils import import_signor, import_kegg, read_uniprot_genename, get_complexes_pairs, get_stringdb, get_stringdb_actions


# -- Imports
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)

# -- Overlap
samples = set(proteomics).intersection(transcriptomics)
proteins = set(proteomics.index).intersection(transcriptomics.index)

# -- Protein-pairs correlation
p_corr = proteomics.ix[proteins, samples].T.corr(method='pearson')
p_corr.values[np.tril_indices(p_corr.shape[0], 0)] = np.nan
p_corr = p_corr.unstack().dropna()

t_corr = transcriptomics.ix[proteins, samples].T.corr(method='pearson')
t_corr.values[np.tril_indices(t_corr.shape[0], 0)] = np.nan
t_corr = t_corr.unstack().dropna()


# -- Interactions
uniprot = read_uniprot_genename()

# CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}

# String
string_thres = {'low': 150, 'medium': 400, 'high': 700, 'highest': 900}

string = {thres: {(uniprot[s][0], uniprot[t][0]) for p1, p2 in get_stringdb(string_thres[thres]) for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot} for thres in string_thres}

string_action = {thres: get_stringdb_actions(string_thres[thres]).groupby('mode')['interaction'].agg(lambda x: set(x)).to_dict() for thres in string_thres}
string_action = {thres: {mode: {(uniprot[s][0], uniprot[t][0]) for p1, p2 in string_action[thres][mode] for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot} for mode in string_action[thres]} for thres in string_thres}

# Signor
signor = import_signor()
signor = {(s, t) for p1, p2 in signor for s, t in [(p1, p2), (p2, p1)]}

# KEGG metabolism
kegg = import_kegg()
kegg = {(s, t) for p1, p2 in kegg for s, t in [(p1, p2), (p2, p1)]}


# -- Plot: Boxplots
labels = {
    'CORUM': 'Complex', 'STRING': 'Functional', 'All': 'All'
}

# All
plot_df = concat({'Proteomics': p_corr, 'Transcriptomics': t_corr}).reset_index()
plot_df.columns = ['Data', 'Px', 'Py', 'Correlation']
plot_df['Interaction'] = 'All'

for n, sets in [('CORUM', corum), ('STRING', string['highest'])]:
    df = concat({'Proteomics': p_corr[sets].dropna(), 'Transcriptomics': t_corr[sets].dropna()}).reset_index()
    df.columns = ['Data', 'Px', 'Py', 'Correlation']
    df['Interaction'] = n

    plot_df = plot_df.append(df, ignore_index=True)

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.violinplot(x='Correlation', y='Interaction', hue='Data', data=plot_df, palette=palette, linewidth=1., cut=0, split=True, inner='quartile', orient='h')
plt.axvline(0.0, ls='--', lw=0.3, c='black', alpha=.5)
plt.xlim(-1, 1)
sns.despine(trim=True)
plt.legend(loc=2)
plt.ylabel('')
g.set_yticklabels([labels[i.get_text()] for i in g.get_yticklabels()])
plt.gcf().set_size_inches(4, .5 * len(set(plot_df['Interaction'])))
plt.savefig('./reports/proteins_correlation_boxplot_all.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/proteins_correlation_boxplot_all.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Correlation distributions: ', './reports/proteins_correlation_boxplot_all.pdf'


# CORUM + STRING
dsets = {
    'CORUM_STRING': [('CORUM', corum), ('STRING', string['highest'])],
    'paths': [('Complex', corum), ('Functional', string['highest']), ('Signalling', signor), ('Metabolism', kegg)],
}

labels = {
    'CORUM_STRING': {'CORUM': 'Direct', 'STRING': 'Functional'}
}

for f in dsets:
    plot_df = DataFrame()

    for n, sets in dsets[f]:
        df = concat({'Proteomics': p_corr[sets].dropna(), 'Transcriptomics': t_corr[sets].dropna()}).reset_index()
        df.columns = ['Data', 'Px', 'Py', 'Correlation']
        df['Interaction'] = n

        plot_df = plot_df.append(df, ignore_index=True)
        print n

    print plot_df

    sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
    g = sns.violinplot(x='Correlation', y='Interaction', hue='Data', data=plot_df, palette=palette, linewidth=1., cut=0, split=True, inner='quartile', orient='h')
    plt.axvline(0.0, ls='--', lw=0.3, c='black', alpha=.5)
    plt.xlim(-1, 1)
    sns.despine(trim=True)
    plt.legend(loc=2)
    g.set_yticklabels([labels[f][i.get_text()] if f in labels else i.get_text() for i in g.get_yticklabels()])
    plt.ylabel('')
    plt.gcf().set_size_inches(4, .5 * len(dsets[f]))
    plt.savefig('./reports/proteins_correlation_boxplot_%s.pdf' % f, bbox_inches='tight')
    plt.savefig('./reports/proteins_correlation_boxplot_%s.png' % f, bbox_inches='tight', dpi=300)
    plt.close('all')
    print '[INFO] Correlation distributions (%s): ' % f, './reports/proteins_correlation_boxplot_%s.pdf' % f

# ROC curves
dsets = {
    'CORUM_STRING': [('CORUM', corum), ('STRING', string['highest'])],
    'paths': [('Complex', corum), ('Functional', string['highest']), ('Signalling', signor), ('Metabolism', kegg)],
}

labels = {
    'CORUM_STRING': {'CORUM': 'Direct', 'STRING': 'Functional'}
}

for f in dsets:
    plot_df = []
    for name, interactions in dsets[f]:
        for dtype, data in [('Proteomics', p_corr), ('Transcriptomics', t_corr)]:
            for i in range(5):
                df = data[interactions].dropna()
                df = DataFrame({'cor': df.append(data.loc[~data.index.isin(interactions)].sample(len(df)))})
                df['tp'] = df.index.isin(interactions).astype(int)

                if len(df['tp']) != 0:
                    curve_fpr, curve_tpr, _ = roc_curve(df['tp'], df['cor'])
                    curve_auc = auc(curve_fpr, curve_tpr)

                    plot_df.append({'Interaction': name, 'Data': dtype, 'AUC': curve_auc, 'FPR': curve_fpr, 'TPR': curve_tpr})

    plot_df = DataFrame(plot_df)

    # Barplot
    sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
    g = sns.factorplot(x='Interaction', y='AUC', data=plot_df, hue='Data', palette=palette, kind='bar', errwidth=.5, legend_out=False)
    plt.axhline(0.5, ls='--', lw=0.3, c='black', alpha=.5)
    # g.despine(left=True)
    g.set_ylabels('AUC')
    plt.ylim(0, 1)
    plt.gcf().set_size_inches(3, 4)
    plt.savefig('./reports/proteins_correlation_roc_%s_barplot.pdf' % f, bbox_inches='tight')
    plt.savefig('./reports/proteins_correlation_roc_%s_barplot.png' % f, bbox_inches='tight', dpi=300)
    plt.close('all')
    print '[INFO] Correlation distributions (%s): ' % f, './reports/proteins_correlation_roc_%s_barplot.pdf' % f

    # ROC curves
    sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
    gs, pos = GridSpec(1, len(dsets[f]), hspace=.3, wspace=.3), 0

    for interaction, _ in dsets[f]:
        ax = plt.subplot(gs[pos])

        df = plot_df[plot_df['Interaction'] == interaction]

        for fpr, tpr, data in df[['FPR', 'TPR', 'Data']].values:
            ax.plot(fpr, tpr, c=sns.light_palette(palette[data], reverse=True).as_hex()[2], lw=.3)

        for dtype in set(df['Data']):
            df_d = df[df['Data'] == dtype]
            df_d = df_d[df_d['AUC'] == df_d['AUC'].median()]

            curve_fpr, curve_tpr, curve_auc = df_d[['FPR', 'TPR', 'AUC']].values[0]

            ax.plot(
                curve_fpr, curve_tpr, c=palette[dtype], lw=1., label='%s (AROC %.2f)' % (dtype, curve_auc),
                path_effects=[pe.Stroke(linewidth=1.2, foreground='white'), pe.Normal()]
            )

        ax.plot([0, 1], [0, 1], 'k--', lw=.3)
        sns.despine(trim=True, ax=ax)

        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.set_title(labels[f][interaction] if interaction in labels else interaction)

        ax.legend(loc='lower right')

        pos += 1

    plt.gcf().set_size_inches(3 * len(dsets[f]), 3)
    plt.savefig('./reports/proteins_correlation_roc_%s.pdf' % f, bbox_inches='tight')
    plt.savefig('./reports/proteins_correlation_roc_%s.png' % f, bbox_inches='tight', dpi=300)
    plt.close('all')
    print '[INFO] Correlation distributions (%s): ' % f, './reports/proteins_correlation_roc_%s.pdf' % f
