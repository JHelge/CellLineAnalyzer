# Load CSV using Pandas
import csv
import pandas
from sklearn.ensemble import RandomForestRegressor
from boruta import BorutaPy


def convert_cnvTable():
    # read in X
    filename = 'WES_pureCN_CNV_genes_20220623.csv'
    X = pandas.read_csv(filename, usecols=["model_name", "symbol", "gene_mean"])
    # remove duplicates
    X = X.drop_duplicates(subset=["model_name", "symbol"], keep='first')
    # cell lines as rows, genes as columns
    X = X.pivot(index='model_name', columns='symbol', values='gene_mean')
    # fill NULL values
    X.fillna(X.mean(), inplace=True)
    # write into file
    X.to_csv('./cnv.csv')
    return X

def convert_results():
    y = pandas.read_excel('GDSC1_fitted_dose_response_24Jul22.xlsx', index_col=0, usecols=["CELL_LINE_NAME", "AUC", "DRUG_ID"])
    y = y.loc[y['DRUG_ID'] == 1047]
    y.drop(['DRUG_ID'], axis=1, inplace=True)
    y.to_csv('auc_results.csv')
    return y

#convert_results()
# COPY NUMBER
#convert_cnvTable()


def convert_expressionTable():
    # read in X
    filename = 'rnaseq_tpm_20220624.csv'
    X = pandas.read_csv(filename)
    print(X.shape)
    X = X.transpose()
    # X.drop(["Unnamed: 1"], axis=0, inplace=True)
    X.columns = X.iloc[0]
    X = X.drop_duplicates(subset=["model_name"])
    X.columns = X.iloc[1]
    X.drop(["model_id", "Unnamed: 1"], axis=0, inplace=True)
    X.set_axis(X.iloc[:, 0], axis=0, inplace=True)
    # X.drop(["model_name"], axis=0, inplace=True)
    X.drop(["symbol"], axis=1, inplace=True)
    X = X.loc[:, X.columns.notna()]
    print(X.shape)
    X.to_csv('exp.csv')
    return X


#convert_expressionTable()
def convert_mutTable():
    X = pandas.read_csv('mutations_summary_20221018.csv', usecols=["gene_symbol", "coding", "model_name"])
    X.drop_duplicates(inplace=True)
   # X = X.pivot(index='model_name', columns='gene_symbol', values='cancer_driver')
    print(X.shape)

    p = X.duplicated(subset=["gene_symbol", "model_name"])
    p = p.to_frame(name="is_dup")
    print((p.loc[p['is_dup'] == True]).shape)
    pX = pandas.merge(p, X, left_index=True, right_index=True)
    pX = pX.loc[pX['is_dup'] == True]
    print(pX.shape)
    pX['coding'] = True
    print(pX.shape)
    X.drop_duplicates(subset=["gene_symbol", "model_name"], keep=False, inplace=True)
    pX = pandas.merge(pX, X, how='outer')
    print(pX.shape)
    X = pX[["gene_symbol", "coding", "model_name"]]

    X["coding"] = X["coding"].astype(float)

    X = X.pivot(index='model_name', columns='gene_symbol', values='coding')

    X.fillna(float(0), inplace=True)
    
    X.to_csv('mut.csv')
    return X


#MUTATION
convert_mutTable()

#FORMAT TABLES
#make tables have same cell lines

X = pandas.read_csv('cnv.csv', index_col=0).add_suffix("_cnv")
print(X.shape)
X.to_csv('cnv_suffix.csv', index=False)



Y = pandas.read_csv('exp.csv', index_col=0).add_suffix("_exp")
print(Y.shape)
Y.to_csv('cnv_rows.csv', index=False)

Z = pandas.read_csv('mut.csv', index_col=0).add_suffix("_mut")
print(Z.shape)
Z.to_csv('mut_rows.csv', index=False)


res = pandas.read_csv('auc_results.csv', index_col=0)
XY = pandas.merge(X, Y, left_index=True, right_index=True)
XYZ = pandas.merge(XY, Z, left_index=True, right_index=True)
XYZr = pandas.merge(XYZ, res, left_index=True, right_index=True)

XYZr.to_csv('merged.csv', index=False)


print(XYZr.shape)

#make files for test run
X = XYZr.filter(like='_cnv').astype(float)
X.dropna(axis=1, how='all')  # drop columns with all values missing
X.fillna(X.mean(), inplace=True)  # fill missing values with column mean
Y = XYZr.filter(like='_exp').astype(float)
Y.dropna(axis=1, how='all')  # drop columns with all values missing
Y.fillna(Y.mean(), inplace=True)  # fill missing values with column mean
Z = XYZr.filter(like='_mut').astype(float)
Z.dropna(axis=1, how='all')  # drop columns with all values missing
Z.fillna(Z.mean(), inplace=True)  # fill missing values with column mean

res = XYZr.filter(items=["AUC"])
res.fillna(res.mean(), inplace=True)  # fill missing values with column mean

X.rename(columns=lambda s: s.replace('_cnv', ''), inplace=True)
Y.rename(columns=lambda s: s.replace('_exp', ''), inplace=True)
Z.rename(columns=lambda s: s.replace('_mut', ''), inplace=True)

print(X.shape)
print(Y.shape)
print(Z.shape)
print(res.shape)

X.to_csv('cnv_rows.csv', index=False)
Y.to_csv('tmp_rows.csv', index=False)
Z.to_csv('mut_rows.csv', index=False)
res.to_csv('res.csv')

#FEATURE SELECTION
# choose subset of samples
#X_new = X.sample(frac = 0.01, axis=1)

#print(X_new.shape)

#X_new = X_new.values
#y = Y.values.ravel()



## define random forest classifier, with utilising all cores and
# sampling in proportion to y labels
#rf = RandomForestRegressor(n_jobs=-1, max_depth=5)

# define Boruta feature selection method
#feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, random_state=1)

# find all relevant features - 5 features should be selected
#feat_selector.fit(X_new, y)

# check selected features - first 5 features are selected
#feat_selector.support_

# check ranking of features
#feat_selector.ranking_

# call transform() on X to filter it down to selected features
#X_filtered = feat_selector.transform(X_new)

#print(X_new)
