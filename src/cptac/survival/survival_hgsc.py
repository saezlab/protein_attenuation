import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import palette_cnv_number, palette
from cptac.utils import gkn, randomise_matrix
from pandas import DataFrame, Series, read_csv, concat
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index, k_fold_cross_validation
from lifelines.statistics import logrank_test
from sklearn.cross_validation import ShuffleSplit, KFold
from sklearn.linear_model.base import LinearRegression
from statsmodels.duration.hazard_regression import PHReg
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict


# -- Clinical
clinical = read_csv('./data/hgsc_clinical.csv', index_col=0)
clinical.index = ['%s-01' % i for i in clinical.index]

clinical['time'] = [i for i in clinical['daystodeath_or_LFU']]
clinical['status'] = [0 if i == 'LIVING' else 1 for i in clinical['vital_status']]

age = read_csv('./data/tcga_clinical.csv', index_col=0)['patient.days_to_birth'].dropna().abs().astype(int)
clinical = clinical[[i[:-3] in age.index for i in clinical.index]]
clinical['age'] = [age.ix[i[:-3]] for i in clinical.index]

clinical = clinical[clinical['time'] < (10 * 365)]
print 'clinical', clinical.shape


# -- Proteomics
proteomics = read_csv('./data/hgsc_proteomics_processed.csv', index_col=0)
proteomics.columns = [i[:15] for i in proteomics]

# Average replicates
remove_samples = {i for i in set(proteomics) if proteomics.loc[:, [i]].shape[1] == 2 and proteomics.loc[:, [i]].corr().ix[0, 1] < .4}
proteomics = proteomics.drop(remove_samples, axis=1)
proteomics = DataFrame({i: proteomics.loc[:, [i]].mean(1) for i in set(proteomics)})

# Drop missing values
proteomics = proteomics.dropna()

# Overlap
proteins, samples = list(set(proteomics.index)), list(set(proteomics).intersection(clinical.index))

# Regress-out age, p = 'EIF3E'
def rm_batch(p):
    ys = proteomics.ix[p, samples].dropna()
    xs = clinical.ix[ys.index, ['age']]

    lm = LinearRegression().fit(xs, ys)

    return ys - xs.dot(lm.coef_) - lm.intercept_

proteomics = DataFrame({p: rm_batch(p) for p in proteins}).T

# Normalise
proteomics = DataFrame({i: gkn(proteomics.ix[i].dropna()).to_dict() for i in proteins}).T

# Export
proteomics.to_csv('./data/hgsc_proteomics_processed_normalised_no_nan.csv')
# proteomics = read_csv('./data/hgsc_proteomics_processed_normalised_no_nan.csv', index_col=0)
print 'proteomics', proteomics.shape


# -- Protein complexes
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if len(corum[k]) > 1}
corum_n = get_complexes_name()
print 'corum', len(corum)


# -- Signif Px -> Py associations
ppairs = read_csv('./tables/ppairs_cnv_regulation.csv')
p_complexes = list({(px, c) for px, py in ppairs[['px', 'py']].values for c in corum if px in corum[c] and py in corum[c]})
print 'p_complexes', len(p_complexes)


# -- Regression: c = 811
alpha, n_folds, n_randomisations = 1., 3, 100

# K-fold
regressions = {c[1]: {} for c in p_complexes}
for permutation in range(15):
    # Complex proteins
    c_proteins = list(corum[c])

    #
    x = concat([proteomics.ix[c_proteins].T, clinical[['time', 'status']]], axis=1).dropna()

    # KFold
    groups, cv = [], KFold(len(samples), n_folds=n_folds, shuffle=True)

    #
    for train, test in cv:
        # Assemble
        x_train, x_test = x.ix[train], x.ix[test]

        # Train: LifeLines
        model = CoxPHFitter(normalize=False, penalizer=alpha).fit(x_train, 'time', event_col='status')

        # Predict
        hazards = model.predict_partial_hazard(x_test[c_proteins])[0]

        # Classify
        groups.append(hazards.apply(lambda i: 1 if i > np.percentile(hazards, 50) else -1))

    # Aggregate predicted samples groups
    groups = concat(groups)

    # Logrank test
    logrank = logrank_test(
        x.ix[groups[groups == 1].index, 'time'], x.ix[groups[groups == -1].index, 'time'],
        x.ix[groups[groups == 1].index, 'status'], x.ix[groups[groups == -1].index, 'status']
    ).test_statistic

    #
    rand_logranks = []
    for randomisation in range(n_randomisations):
        rand_x = concat([proteomics.ix[c_proteins].T, randomise_matrix(clinical['time']), randomise_matrix(clinical['status'])], axis=1).dropna()

        rand_groups = []
        for train, test in cv:
            # Assemble
            x_train, x_test = rand_x.ix[train], rand_x.ix[test]

            # Train: LifeLines
            model = CoxPHFitter(normalize=False, penalizer=alpha).fit(x_train, 'time', event_col='status')

            # Predict
            hazards = model.predict_partial_hazard(x_test[c_proteins])[0]

            # Classify
            rand_groups.append(hazards.apply(lambda i: 1 if i > np.percentile(hazards, 50) else -1))

        # Aggregate predicted samples groups
        rand_groups = concat(rand_groups)

        # Logrank test
        rand_logranks.append(logrank_test(
            rand_x.ix[rand_groups[rand_groups == 1].index, 'time'], rand_x.ix[rand_groups[rand_groups == -1].index, 'time'],
            rand_x.ix[rand_groups[rand_groups == 1].index, 'status'], rand_x.ix[rand_groups[rand_groups == -1].index, 'status']
        ).test_statistic)

    p_value = float(sum([logrank > i for i in rand_logranks])) / n_randomisations
    p_value = 1. / n_randomisations if p_value == 0 else p_value

    print logrank, p_value

    # Store results
    regressions[c][permutation] = groups


# -- Plot
c = 811

df = DataFrame(regressions[c]).mode(axis=1).apply(lambda x: min(x.min(), x.max(), key=abs), axis=1)
df = concat([df, clinical[['time', 'status']]], axis=1).dropna()

samples_up = set(df[df[0] == 1].index)
samples_down = set(df[df[0] == -1].index)

logrank = logrank_test(
    df.ix[samples_up, 'time'], df.ix[samples_down, 'time'],
    df.ix[samples_up, 'status'], df.ix[samples_down, 'status']
)

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'}, font_scale=0.75)
ax = plt.subplot(111)

kmf = KaplanMeierFitter()

kmf.fit(df.ix[samples_up, 'time'], event_observed=df.ix[samples_up, 'status'], label='Up')
kmf.plot(ci_force_lines=False, color=palette_cnv_number[1], ax=ax, show_censors=True, censor_styles={'ms': 5})

kmf.fit(df.ix[samples_down, 'time'], event_observed=df.ix[samples_down, 'status'], label='Down')
kmf.plot(ci_force_lines=False, color=palette_cnv_number[-1], ax=ax, show_censors=True, censor_styles={'ms': 5})

sns.despine(ax=ax)

ax.set_title('%s\nLogrank test: %.2e' % (corum_n[c], logrank.p_value))
ax.set_xlabel('Timeline (days)')
ax.set_ylabel('Survival fraction')

plt.gcf().set_size_inches(5, 3)
plt.savefig('./reports/survival_protein_complex.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
