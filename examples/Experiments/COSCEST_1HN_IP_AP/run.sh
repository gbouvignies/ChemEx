#!/bin/bash

SELECTION="48 61 114 115 147 151"

chemex fit -e Experiments/*.toml \
           -p Parameters/*.toml \
           -m Methods/*.toml \
           -d 2st_rs \
           --include $SELECTION \
           -o Output
