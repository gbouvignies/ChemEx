#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/global.toml \
              Parameters/n15_cs.toml \
              Parameters/n15_dw.toml \
           -m Methods/method.toml \
           -d 2st.pb_kex \
           -o Output
