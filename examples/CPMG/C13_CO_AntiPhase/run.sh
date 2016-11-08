#!/bin/sh

chemex fit -e Experiments/13co_ap_*.cfg \
           -p Parameters/params_co.cfg \
           -m Methods/method.cfg \
           -d 2st.pb_kex \
           -o Output
