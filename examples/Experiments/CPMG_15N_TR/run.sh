#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/global.toml \
              Parameters/cs_a.toml \
              Parameters/dw_ab.toml \
              Parameters/r1*.toml \
           -m Methods/method.toml \
           -o Output
