![Python tests](https://github.com/jodyphelan/tb-ml/workflows/Tests/badge.svg)

# tb-ml

A standardised approach to predicting antimicrobial resistance from WGS data using machine learning models or direct association. The pipeline essentially runs two user-defined docker containers: one for data variant calling (VC) and one to perform resistance prediction. We provide example containers for both steps (<https://github.com/julibeg/tb-ml-containers> and <https://hub.docker.com/u/julibeg>). The idea is to allow users to create / swap out different containers for both steps and thus compare the results of different prediction models (or VC pipelines).

New containers created by researchers developing new prediction models should adhere to the following standards:

* The **prediction container** should be able to be called with one of two arguments: `get_target_vars` or `predict`. In the first case (i.e. called as `docker run image-name get_target_vars`), the container should print a CSV with the target variants it expects for prediction in the format `POS,REF,ALT,AF` (`AF` being the allele frequency of the training dataset) to `STDOUT`. In the second case, (i.e. called as `docker run image-name predict input.csv`), it should read the input file (in the format `POS,REF,ALT,GT`) and print the prediction result (a single floating point number between `0` and `1`) to `STDOUT`.

* The **VC container** should be called like this: `docker run image-name -b path-to-BAM-file -t path-to-CSV-with-target-variants`. It should accept an alignment file (SAM/BAM/CRAM) and a CSV in the format `POS,REF,ALT,AF` as inputs. The purpose of the CSV is to make sure that the positions required by the prediction container are covered in the called variants.

## Install

```bash
git clone https://github.com/jodyphelan/tb-ml.git
pip install -r requirements.txt 
pip install .
```
