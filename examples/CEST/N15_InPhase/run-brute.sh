#!/bin/sh

chemex fit -e Experiments/cest_n15*.cfg \
           -p Parameters/params_n15_brute.cfg \
           -m Methods/method_n15_brute_indiv.cfg \
           -d 2st.pb_kex \
           +r 50N-HN \
           -o Output/brute

chemex fit -e Experiments/cest_n15*.cfg \
           -p Parameters/params_n15_brute.cfg \
           -m Methods/method_n15_brute_global.cfg \
           -d 2st.pb_kex \
           -o Output/brute-global
