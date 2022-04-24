#!/bin/bash

bam_file=$PWD/../../test_data/test.cram

docker build -t vc-test --rm . > /dev/null

bash dev-get-target-vars.sh | docker run -i \
    --mount type=bind,source="$bam_file",target=/data/aligned_reads \
    vc-test
