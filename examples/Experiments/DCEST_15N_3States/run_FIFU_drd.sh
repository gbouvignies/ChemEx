#!/bin/sh

# Correct FIFU model, including DRD-CEST datasets
chemex fit -e Experiments/*.toml \
              Experiments/DRD/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -d 3st_fork \
           -o Output_FIFU_DRD
