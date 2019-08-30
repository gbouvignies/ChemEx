#!/bin/sh

chemex fit -e Experiments/20Hz_?00Hz.toml \
           -p Parameters/global.toml \
              Parameters/cs_a.toml \
              Parameters/dw_ab.toml \
           -m Methods/method.toml \
           -o Output
