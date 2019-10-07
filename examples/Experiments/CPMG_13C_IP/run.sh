#!/bin/sh

chemex fit \
    -e Experiments/*.toml \
    -p Parameters/global.toml \
       Parameters/cs_a.toml \
    -m Methods/method.toml \
    -o Output
