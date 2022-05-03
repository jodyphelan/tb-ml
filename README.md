![Python tests](https://github.com/jodyphelan/tb-ml/workflows/Tests/badge.svg)

# tb-ml

A standardised approach to predicting antimicrobial resistance from WGS data using machine learning models or direct association. The pipeline essentially runs two user-defined docker containers: one for variant calling and one to perform resistance prediction. We provide default containers for both steps on Docker Hub (at https://hub.docker.com/u/julibeg and https://hub.docker.com/u/jodyphelan; the containers were built with files in the `/docker` directory of this repo). This allows users to swap different containers for both steps and compare the results of different models if containers are available for them. 

For now, the variant calling container is expected to read a file with aligned reads (`.sam`, `.sam`, or `.cram`) from `/data/aligned_reads` and a CSV with the variants of interest in the format `POS,REF,ALT` (**with** a header line) from `STDIN` and print the results as `POS,REF,ALT,GT,DP` (**with** a header line) to `STDOUT`. We will add support for containers that take raw reads and run the complete variant calling pipeline in the near future. Non-calls will be replaced with the AF in the training dataset by `tb-ml` before calling the prediction container.

The prediction container is expected to read a CSV in the form of `POS,REF,ALT,GT` (**with** a header line) from `STDIN` and print the prediction result for the sample to `STDOUT`. `tb-ml` will read the result and print a short comma-separated report. The prediction container also needs to return the allele frequencies of the dataset it has been trained on in the format `POS,REF,ALT,AF` (**with** a header line) when called with `get_target_vars` as single argument. 

## Install

```bash
git clone https://github.com/jodyphelan/tb-ml.git
pip install .
```
