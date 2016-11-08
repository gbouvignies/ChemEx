#!/bin/sh

chemex fit -e Experiments/h1_ap_*.cfg \
           -p Parameters/params_h1.cfg \
           -m Methods/method_h1.cfg \
           -d 2st.pb_kex \
           -o Output
