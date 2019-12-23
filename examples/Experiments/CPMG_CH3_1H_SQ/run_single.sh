#!/bin/sh

# Optional, perform single-residue fit first
chemex fit -e Experiments/*.toml \
           -p Parameters/*.toml \
           -d 2st_rs \
           -o Output/Single
