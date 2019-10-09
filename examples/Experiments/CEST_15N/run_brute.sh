#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/global_brute.toml \
              Parameters/cs_a.toml \
              Parameters/dw_ab_brute.toml \
           -m Methods/method_brute_indiv.toml \
           -o Output/Brute \
           -r 50N

chemex fit -e Experiments/*.toml \
           -p Parameters/global_brute.toml \
              Parameters/cs_a.toml \
              Parameters/dw_ab_brute.toml \
           -m Methods/method_brute_global.toml \
           -o Output/BruteGlobal
