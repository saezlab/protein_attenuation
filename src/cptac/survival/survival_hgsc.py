import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import palette_cnv_number, palette
from cptac.utils import gkn
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
p_complexes = {(px, c) for px, py in ppairs[['px', 'py']].values for c in corum if px in corum[c] and py in corum[c]}
print 'p_complexes', len(p_complexes)


# --
# px, c = 'EIF3A', 1097
c_survival = {}
for px, c in p_complexes:
    c_proteins = list(corum[c])
    print corum_n[c], len(c_proteins)

    # Assemble data-set
    df = concat([proteomics.ix[c_proteins].T, clinical[['time', 'status']]], axis=1).dropna()

    # Bootstrap
    # cv = ShuffleSplit(len(df), test_size=.1, n_iter=20)
    cv = KFold(len(df), n_folds=5)

    cindexes_train, cindexes_test, group = [], [], []
    for train, test in cv:
        x_train, x_test = df.ix[train], df.ix[test]

        # # -- Statsmodels
        # # Train
        # model = PHReg(x_train['time'], x_train[c_proteins], x_train['status'], ties='efron')
        # model = model.fit_regularized(L1_wt=0, alpha=.1)
        # # print model.summary()
        #
        # # Predict
        # hazards_train = Series(dict(zip(*(x_train.index, model.predict(exog=x_train[c_proteins], pred_type='hr').predicted_values))))
        # hazards_test = Series(dict(zip(*(x_test.index, model.predict(exog=x_test[c_proteins], pred_type='hr').predicted_values))))

        # -- LifeLines
        # Train
        model = CoxPHFitter(normalize=False, penalizer=1.)
        model = model.fit(x_train, 'time', event_col='status')
        # print model.print_summary()

        # Predict
        hazards_train = model.predict_partial_hazard(x_train[c_proteins])[0]
        hazards_test = model.predict_partial_hazard(x_test[c_proteins])[0]

        # -- Evalute
        # Concordance index (c-index)
        cindex_train = concordance_index(x_train['time'], -hazards_train, x_train['status'])
        cindex_test = concordance_index(x_test['time'], -hazards_test, x_test['status'])

        # Classify into groups
        groups_train = hazards_train.apply(lambda i: 1 if i > np.percentile(hazards_train, 50) else -1)
        groups_test = hazards_test.apply(lambda i: 1 if i > np.percentile(hazards_test, 50) else -1)

        # # Log-rank test
        # rank_train = logrank_test(
        #     x_train.ix[groups_train[groups_train == 1].index, 'time'], x_train.ix[groups_train[groups_train == -1].index, 'time'],
        #     x_train.ix[groups_train[groups_train == 1].index, 'status'], x_train.ix[groups_train[groups_train == -1].index, 'status']
        # )
        #
        # rank_test = logrank_test(
        #     x_test.ix[groups_test[groups_test == 1].index, 'time'], x_test.ix[groups_test[groups_test == -1].index, 'time'],
        #     x_test.ix[groups_test[groups_test == 1].index, 'status'], x_test.ix[groups_test[groups_test == -1].index, 'status']
        # )

        # -- Store results
        group.append(groups_test)
        cindexes_train.append(cindex_train)
        cindexes_test.append(cindex_test)

        # print 'cindex: train %.2f, test %.2f' % (cindex_train, cindex_test)
        # print 'logrank: train %.2e, test %.2f' % (rank_train.p_value, rank_test.p_value)

    cindexes_train_mean, cindexes_test_mean = np.median(cindexes_train), np.median(cindexes_test)
    print 'train %.2f, test %.2f\n' % (cindexes_train_mean, cindexes_test_mean)

    c_survival[c] = {
        'complex': c, 'len': len(c_proteins), 'group': concat(group),
        'cindexes_train': cindexes_train, 'cindexes_test': cindexes_test,
        'cindexes_train_mean': cindexes_train_mean, 'cindexes_test_mean': cindexes_test_mean
    }

c_survival = DataFrame(c_survival).T
c_survival['name'] = [corum_n[i] for i in c_survival.index]
print c_survival.sort('cindexes_test_mean', ascending=False)[['cindexes_train_mean', 'cindexes_test_mean', 'name', 'len']]


# -- Plot
c = 811

df = c_survival.ix[c, 'group']
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
