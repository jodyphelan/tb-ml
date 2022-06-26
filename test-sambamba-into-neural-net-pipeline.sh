#!/bin/bash

proc_cont="julibeg/tb-ml-one-hot-encoded-from-cram:v0.3.0"
pred_cont="julibeg/tb-ml-neural-net-predictor-13-drugs:v0.4.0"

tb-ml \
    --container $pred_cont \
        "--get-target-loci \
        -o target-loci.csv" \
    --container $proc_cont \
        "-b tb_ml/tests/test_data/test.cram \
        -r target-loci.csv \
        -o one-hot-seqs.csv" \
    --container $pred_cont \
        "one-hot-seqs.csv"