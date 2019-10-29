#!/bin/sh

chemex simulate \
    -e Experiments/*.toml \
    -p Parameters/*.toml \
    -o OutputSim
