![Python tests](https://github.com/jodyphelan/tb-ml/workflows/Tests/badge.svg)

# tb-ml

A plethora of machine learning models for resistance prediction from WGS data are being published every year, but in many cases neither the source code nor the final model are provided, making systematic comparisons difficult. This package introduces a standardised antimicrobial resistance prediction pipeline using Docker containers with either machine learning models or methods for direct association. The pipeline essentially runs two Docker containers specified by the user: one for variant calling (VC) and one to perform resistance prediction (example containers can be found on [Docker Hub](https://hub.docker.com/u/julibeg) and [GitHub](https://github.com/julibeg/tb-ml-containers)). We define a flexible standard for the call signatures of both containers (for details see below). As long as they adhere to the standard, users can create and swap out containers for both steps, effectively facilitating easy comparisons between different methods for resistance prediction.

## Install

```bash
git clone https://github.com/jodyphelan/tb-ml.git
pip install -r requirements.txt 
pip install .
```

## Example usage

```bash
tb-ml -b tb_ml/tests/test_data/test.cram \
    -v julibeg/tb-ml-freebayes-vc-from-cram:v0.1.0 \
    -p julibeg/tb-ml-simple-rf-predictor-streptomycin:v0.1.0
```

This command specifies the input CRAM file and the names of the containers to use for the variant calling (`-v`) and prediction (`-p`) steps. It runs the prediction and prints a short report to `STDOUT`:

```
parameter,value
file,/abs/path/to/tb-ml/tb_ml/tests/test_data/test.cram
vc_container,julibeg/tb-ml-freebayes-vc-from-cram:v0.1.0
pred_container,julibeg/tb-ml-simple-rf-predictor-streptomycin:v0.1.0
shared_variants,7248
dropped_variants,0
missing_variants,28
noncalls,96
variants_set_to_AF,124
resistance_probability,0.33
resistance_status,S
```

## New containers

Containers created by researchers developing new prediction models should adhere to the standards outlined below in order to be used in `tb-ml`. All the pre-processing should take place in the VC container and the prediction container should only take a vector of features and produce a single number: the probability of resistance. One peculiarity of the prediction pipeline is that the VC container needs to provide variant calls for all target variants used by the prediction container. Therefore, the prediction container must implement a method for giving out a list of the variants which can be used by the VC container to make sure that all relevant variants are present in its output. The call signatures are as follows:

* The **prediction container** should be called with one of two arguments: `get_target_vars` or `predict`. In the first case (i.e. called as `docker run image-name get_target_vars`), the container should print a CSV with the target variants it expects for prediction in the format `POS,REF,ALT,AF` (`AF` being the allele frequency of the training dataset) to `STDOUT`. In the second case, (i.e. called as `docker run image-name predict input.csv`), it should read the input file (in the format `POS,REF,ALT,GT`) and print the prediction result (a single floating point number between `0` and `1`) to `STDOUT`.

* The **VC container** should be called like this: `docker run image-name -b path-to-BAM-file -t path-to-CSV-with-target-variants`. It should accept an alignment file (SAM/BAM/CRAM) and a CSV in the format `POS,REF,ALT,AF` as inputs. The purpose of the CSV is to make sure that the positions required by the prediction container are covered in the called variants. Non-calls and variants missing from the VC output should be set to the corresponding allele frequencies (available from the input file with the target variants). The output should be written to `STDOUT` and be in the format `POS,REF,ALT,GT`. It can start with comments with some basic info / stats about the variant calling process. The comment lines should start with `#` and hold comma-separated key-value pairs. They will be added to the prediction report.
