#!/bin/sh

chemex fit -e Experiments/*.cfg \
           -p Parameters/params_n15.cfg \
           -d 2st.pb_kex \
           -o Output
