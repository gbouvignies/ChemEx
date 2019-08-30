#!/bin/bash

chemex fit \
    -e Experiments/*.toml \
    -p Parameters/global.toml \
       Parameters/cs_a.toml \
       Parameters/dw_ab.toml \
    -m Methods/method.toml \
    -o Output
