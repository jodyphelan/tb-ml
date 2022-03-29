import joblib
import pandas as pd
from sklearn import metrics as mets

m = joblib.load('RF_INH_fitted.pkl')
xtrain, xtest, ytrain, ytest = joblib.load('INH_train_test_split.pkl')

ypred = m.predict(xtest)
print(mets.f1_score(ytest, ypred))