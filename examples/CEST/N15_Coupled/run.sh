#!/bin/sh

chemex fit -e Experiments/cest_n15_500ms_25Hz.cfg \
           -p Parameters/params_n15.cfg \
           -m Methods/method_n15.cfg \
           -d 2st.pb_kex \
           -o Output
