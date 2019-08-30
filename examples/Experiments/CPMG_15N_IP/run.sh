#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/global.toml \
              Parameters/cs_a.toml \
              Parameters/dw_ab.toml \
              Parameters/r1_a.toml \
           -m Methods/method.toml \
           -d 2st.pb_kex \
           -o Output
