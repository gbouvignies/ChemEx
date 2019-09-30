#!/bin/bash

rm -rf Output/Global

chemex fit \
    -e Experiments/*.toml \
    -p Parameters/*.toml \
    -m Methods/method.toml \
    -o Output/Global
