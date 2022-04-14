#!/bin/bash

bam_file=$PWD/../../test_data/test.cram

docker build -t vc-test --rm . > /dev/null

python dev_get_target_vars.py | docker run -i \
    --mount type=bind,source="$bam_file",target=/data/aligned_reads \
    vc-test
