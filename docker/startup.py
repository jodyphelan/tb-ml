import joblib
import pandas as pd

m = joblib.load('../test_models/RF_INH_fitted.pkl')
print(m.feature_importances_[:10])