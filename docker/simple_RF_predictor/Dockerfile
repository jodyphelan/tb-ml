FROM python:3.8-bullseye

RUN pip install numpy pandas joblib scikit-learn

# copy the python script and the model / data files
COPY script.py /
COPY SM_RF_fitted.pkl /model.pkl
COPY SM_RF_target_vars_AFs.csv /target_vars.csv

ENTRYPOINT ["python", "/script.py"]
