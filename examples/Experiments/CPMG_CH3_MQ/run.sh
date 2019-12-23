#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/*.toml \
           -m Methods/*.toml \
           -o Output
