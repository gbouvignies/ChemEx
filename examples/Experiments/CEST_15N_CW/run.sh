#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/global.toml \
              Parameters/cs_a_*.toml \
              Parameters/dw_ab_*.toml \
              Parameters/r1a_a.toml \
           -m Methods/method.toml \
           -o Output
