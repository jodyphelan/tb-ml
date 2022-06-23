#!/bin/bash

cont1="julibeg/tb-ml-one-hot-encoded-from-cram:v0.2.0"
cont2="julibeg/tb-ml-neural-net-predictor-13-drugs:v0.3.0"

tb-ml \
    --container $cont1 \
    "-b tb_ml/tests/test_data/test.cram \
        -r /home/julian/git/tb-ml-containers/one_hot_encode/dev/H37RV_loci_coords_from_table_in_paper.bed \
        -o one-hot-seqs.csv" \
    --container $cont2 \
        one-hot-seqs.csv