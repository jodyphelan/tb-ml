import setuptools

version = [l.strip() for l in open("tb_ml/__init__.py") if
           "version" in l][0].split('"')[1]

setuptools.setup(

    name="tb-ml",
    version=version,
    packages=["tb_ml"],
    license="GPL3",
    long_description="A standardised approach to predicting resistance using machine learning",
    scripts=[
        'scripts/tb-ml-predict.py'
    ]
)