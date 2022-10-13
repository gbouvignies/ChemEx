#!/bin/sh

# Incorrect FIU model, including DRD-CEST datasets
chemex fit -e Experiments/*.toml \
              Experiments/DRD/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -d 3st_linear \
           -o Output_FIU_DRD
