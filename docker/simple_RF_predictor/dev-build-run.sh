#!/bin/bash

docker build -t rf-sm-test --rm .

docker run -i rf-sm-test \
    predict \
    < ../../test_data/test_variants.csv