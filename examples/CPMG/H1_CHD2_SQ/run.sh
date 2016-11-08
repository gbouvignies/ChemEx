#!/bin/sh

chemex fit -e Experiments/chd2_h1sq_*.cfg \
           -p Parameters/params_chd2_h1sq.cfg \
           -m Methods/method_chd2.cfg \
           -d 2st.pb_kex \
           -o Output
