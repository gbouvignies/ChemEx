#!/bin/sh

chemex fit -e Experiments/*hz_fast.toml \
           -p Parameters/global.toml \
              Parameters/n15_cs.toml \
              Parameters/n15_dw.toml \
           -m Methods/method.toml \
           -o Output/Fast
