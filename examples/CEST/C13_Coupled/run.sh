#!/bin/sh

chemex fit -e Experiments/cest_c13_ali_250ms_23Hz*cfg \
           -p Parameters/params_c13.cfg \
           -m Methods/method_c13.cfg \
           -d 2st.pb_kex \
           -o Output
