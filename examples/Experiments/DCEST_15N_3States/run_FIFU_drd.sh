#!/bin/sh

SELECTION="19 20 26 32 35 39"

# Correct FIFU model, including DRD-CEST datasets
chemex fit -e Experiments_DRD/*.toml \
           -p Parameters_FIFU/*.toml \
           -m Methods_FIFU/*.toml \
           -d 3st \
           --include $SELECTION \
           -o Output_FIFU_DRD
