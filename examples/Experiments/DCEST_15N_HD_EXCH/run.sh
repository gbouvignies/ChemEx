#!/bin/sh

chemex fit \
    -e Experiments/*.toml \
    -p Parameters/global.toml \
       Parameters/cs_a.toml \
    -m Methods/method.toml \
    -d 2st.hd_exch \
    -o Output \
    -r 57N 66N 74N 78N 79N 93N 96N 97N 98N 109N
