import setuptools

with open("tb_ml/__init__.py") as f:
    version = [line.strip()
               for line in f if "version" in line][0].split('"')[1]

setuptools.setup(
    name="tb-ml",
    version=version,
    packages=["tb_ml"],
    license="GPL3",
    long_description="A standardised approach to predicting resistance using' \
         machine learning",
    entry_points={
        "console_scripts": [
            'tb-ml = tb_ml.__main__:main',
        ]
    },
)
