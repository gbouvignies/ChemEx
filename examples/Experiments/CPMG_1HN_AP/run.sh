#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/global.toml \
              Parameters/cs_a.toml \
              Parameters/r1_a*.toml \
              Parameters/r1a_a.toml \
           -m Methods/method.toml \
           -o Output
