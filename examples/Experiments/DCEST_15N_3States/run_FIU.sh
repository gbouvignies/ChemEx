#!/bin/sh

# Incorrect FIU model, without DRD-CEST datasets
chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -d 3st_linear \
           -o Output_FIU
