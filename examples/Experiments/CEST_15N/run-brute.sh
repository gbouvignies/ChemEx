#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/global_brute.toml \
              Parameters/n15_cs.toml \
              Parameters/n15_dw_brute.toml \
           -m Methods/method_brute_indiv.toml \
           -o Output/Brute \
           +r 50N-HN

chemex fit -e Experiments/*.toml \
           -p Parameters/global_brute.toml \
              Parameters/n15_cs.toml \
              Parameters/n15_dw_brute.toml \
           -m Methods/method_brute_global.toml \
           -o Output/BruteGlobal
