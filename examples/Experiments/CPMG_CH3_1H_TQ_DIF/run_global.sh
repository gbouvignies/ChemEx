#!/bin/bash

rm -rf output_global

chemex fit -e Experiments/*.toml \
           -p Parameters/*.toml \
           -m Methods/method_global.toml \
           -o Output
