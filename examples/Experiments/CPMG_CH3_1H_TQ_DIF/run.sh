#!/bin/bash

chemex fit -e Experiments/*.toml \
           -p Parameters/*.toml \
           -m Methods/method.toml \
           -o Output
