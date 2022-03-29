import numpy as np
import pandas as pd
import joblib
import sys

m = joblib.load('model.pkl')

X = pd.read_csv(sys.stdin, index_col=0)

if X.shape[1] == 1:
    # the CSV had only 2 columns --> there is only one sample
    X = X.T

ypred = m.predict(X)
print(np.squeeze(ypred))
