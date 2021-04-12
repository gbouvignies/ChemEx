#!/bin/sh

SELECTION="19 20 26 32 35 39"

# Incorrect FIU model, including DRD-CEST datasets
chemex fit -e Experiments/*.toml \
              Experiments/DRD/*.toml \
           -p Parameters/parameters_FIU.toml \
           -m Methods/method_FIU.toml \
           -d 3st \
           --include $SELECTION \
           -o Output_FIU_DRD
