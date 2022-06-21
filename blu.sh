#!/bin/bash

pred_container="rf-test"
vc_container="vc-test"
# pred_container="julibeg/tb-ml-simple-rf-predictor-streptomycin:v0.2.0"
# vc_container="julibeg/tb-ml-freebayes-vc-from-cram:v0.1.0"

tb-ml \
    --container $pred_container \
        "--get-target-vars \
        -o target_vars.csv" \
    --container $vc_container \
        "-b tb_ml/tests/test_data/test.cram \
        -t target_vars.csv \
        -o genotypes.csv" \
    --container $pred_container \
        "genotypes.csv \
        -o prediction_result.csv"

