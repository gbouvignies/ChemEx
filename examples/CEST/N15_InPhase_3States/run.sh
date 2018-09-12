#!/bin/sh

chemex fit -e experiments/cest_n15*.cfg \
           -p parameters/params_n15.cfg \
           -m methods/method_n15.cfg \
           -d 3st.pb_kex \
           -o Output
