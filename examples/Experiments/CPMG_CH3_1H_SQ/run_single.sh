#!/bin/sh

# Optional, perform single-residue fit first
chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -d 2st.rs \
           -o Output/Single
