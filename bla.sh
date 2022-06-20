#!/bin/bash

tb-ml --preprocessing-container julibeg/tb-ml-freebayes-vc-from-cram:v0.1.0 \
    --preprocessing-args "" \
    --prediction-container julibeg/tb-ml-simple-rf-predictor-streptomycin:v0.2.0
    --prediction-args ""