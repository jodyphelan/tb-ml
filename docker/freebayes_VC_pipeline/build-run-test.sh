#!/bin/bash

bam_file=$PWD/test_data/test.cram
target_vars=$PWD/test_data/target_vars.csv

docker build . -t vc-test --rm

docker run -i \
    --mount type=bind,source="$bam_file",target=/data/aligned_reads \
    vc-test < "$target_vars"
