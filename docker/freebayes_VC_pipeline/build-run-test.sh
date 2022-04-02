#!/bin/bash

docker build . -t vc-test --rm && docker run -it \
    -v /home/julsy/git/tb-ml/docker/freebayes_VC_pipeline/data:/data \
    -v $PWD/bla:/bla \
    vc-test \
    SM_training_AF.csv \
    test.cram 