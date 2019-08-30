#!/bin/bash

chemex fit \
    -e Experiments/*.toml \
    -p parameters/global.toml \
       parameters/cs_a.toml \
    -m methods/method.toml \
    -d 2st.pb_kexrs \
    -o Output \
    +r 11N-HN 25N-HN 60N-HN 62N-HN \
       63N-HN 66N-HN 80N-HN 81N-HN \
       82N-HN 83N-HN 96N-HN
