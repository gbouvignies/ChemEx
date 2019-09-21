#!/bin/bash

chemex fit \
    -e Experiments/*.toml \
    -p Parameters/global.toml \
       Parameters/cs_a.toml \
    -m Methods/method.toml \
    -d 2st.pb_kexrs \
    -o Output \
    -r 11N 25N 60N 62N 63N 66N 80N 81N 82N 83N 96N
