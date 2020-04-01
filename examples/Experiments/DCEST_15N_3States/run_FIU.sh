#!/bin/sh

SELECTION="19 20 26 32 35 39"

# Incorrect FIU model, without DRD-CEST datasets
chemex fit -e Experiments/*.toml \
           -p Parameters/cs_a.toml \
              Parameters/dw_*.toml \
              Parameters/global_FIU.toml \
           -m Methods/method_FIU.toml \
           -d 3st \
           --include $SELECTION \
           -o Output_FIU
