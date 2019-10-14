#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/*.toml \
           -o Output
