FROM continuumio/miniconda3:latest

RUN conda install -c conda-forge -y \
    numpy=1.22.3 \
    pandas=1.4.1 \
    joblib=1.1.0 \
    scikit-learn=1.0.2

# copy the python script and the model / data files
COPY script.py /
COPY SM_RF_fitted.pkl /model.pkl
COPY SM_RF_target_vars_AFs.csv /target_vars.csv

ENTRYPOINT ["python", "/script.py"]
