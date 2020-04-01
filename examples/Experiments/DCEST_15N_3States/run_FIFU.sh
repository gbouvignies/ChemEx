#!/bin/sh

SELECTION="19 20 26 32 35 39"

# Correct FIFU model, without DRD-CEST datasets
chemex fit -e Experiments/*.toml \
           -p Parameters/cs_a.toml \
              Parameters/dw_*.toml \
              Parameters/global_FIFU.toml \
           -m Methods/method_FIFU.toml \
           -d 3st \
           --include $SELECTION \
           -o Output_FIFU
