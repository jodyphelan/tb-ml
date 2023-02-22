![Python tests](https://github.com/jodyphelan/tb-ml/workflows/Tests/badge.svg)

# tb-ml

Large numbers of machine learning models for resistance prediction from WGS data are being published every year, but in many cases neither the source code nor the final model are provided, making systematic comparisons difficult. This package introduces a standardised antimicrobial resistance prediction pipeline using Docker containers with either machine learning models or methods for direct association. The pipeline essentially runs multiple Docker containers specified by the user in succession. Usually, this would include one container for pre-processing of the genomic data (e.g. variant calling) and one to perform resistance prediction (example containers can be found on [Docker Hub](https://hub.docker.com/u/julibeg) and [GitHub](https://github.com/julibeg/tb-ml-containers)). We define a flexible standard for the call signatures of the containers (for details see below). As long as they adhere to the standard, users can create and swap out containers for steps of the pipeline, effectively facilitating easy comparisons between different methods for resistance prediction.

## Install

`tb-ml` requires [Docker](https://www.docker.com/) to be installed and [set up to be run without sudo](https://docs.docker.com/engine/install/linux-postinstall/). If you have Docker, install `tb-ml` with

```bash
git clone https://github.com/jodyphelan/tb-ml.git
cd tb-ml
pip install .
```

## Usage

`tb-ml` can take an arbitrary number of Docker containers and accompanying arguments. It contains a basic Docker API and runs the containers in succession. The CLI only has a single flag (`--container`). Use it to specify the name of a docker image and add a string holding the arguments to be passed to the container after that.

The following example uses a [neural network container](https://github.com/julibeg/tb-ml-containers/tree/main/neural_net_predictor_13_drugs) for prediction, which accepts one-hot-encoded sequences for a number of genomic loci as input. These one-hot-encoded sequences are prepared by a [pre-processing container](https://github.com/julibeg/tb-ml-containers/tree/main/one_hot_encode).In the example, we use a container that extracts the consensus sequence for an Mtb sample from the corresponding SAM/BAM file holding reads aligned against the reference genome H37Rv (ASM19595v2). The pre-processing container expects a list of target loci and the prediction container implements a method to get a list of these loci. Therefore, the pipeline run by `tb-ml` will comprise three steps:

* Get the target loci from the **prediction container**.
* Give the target loci and input alignment file to the **pre-processing container** to
  produce the one-hot-encoded sequences.
* Feed the one-hot-encoded sequences to the **prediction container** to predict the
  resistance status of the sample.

To execute these steps with `tb-ml`, run the following from the root directory of the git repository

```bash
proc_cont="julibeg/tb-ml-one-hot-encoded-seqs-from-raw-reads:v0.2.0"
pred_cont="julibeg/tb-ml-neural-net-from-one-hot-encoded-seqs-13-drugs:v0.7.0"

tb-ml \
    --container $pred_cont \
    "--get-target-loci \
        -o target-loci.csv" \
    --container $proc_cont \
    "-r target-loci.csv \
        -o one-hot-seqs.csv \
        tb_ml/tests/test_data/test_raw_reads_1.fastq.gz \
        tb_ml/tests/test_data/test_raw_reads_2.fastq.gz" \
    --container $pred_cont \
    "one-hot-seqs.csv"
```

## New containers

`tb-ml` was designed with maximum flexibility in mind. New containers only need to adhere to the following rules:

* Any output that should be added to the final report should be printed to `STDOUT`.

* Output that will be used by other containers needs to be written to files.

* Prediction containers should **only predict** and do no pre-processing. In most cases, they should accept input in CSV format (e.g. one-hot-encoded sequences, called variants, etc.)

* Prediction containers can have extra methods to give information that might be relevant for pre-processing (e.g. target loci, variants, etc.). Since this information is used by other containers, it should be written (ideally as CSV) to a file.
