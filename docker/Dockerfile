FROM python:3.8-bullseye

RUN pip install numpy pandas joblib scikit-learn

# copy the python script and the model / data files
COPY startup.py /
COPY *.pkl /

ENTRYPOINT ["python", "./startup.py"]
