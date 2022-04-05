#!/bin/bash

af_file=$PWD/test_data/../test_data/SM_training_AF.csv
bam_fname=test.cram
bam_file=$PWD/test_data/$bam_fname

docker build . -t vc-test --rm

docker run \
    --cpus=6 \
    --mount type=bind,source="$af_file",target=/data/AFs.csv \
    --mount type=bind,source="$bam_file",target=/data/$bam_fname \
    vc-test \
    $bam_fname