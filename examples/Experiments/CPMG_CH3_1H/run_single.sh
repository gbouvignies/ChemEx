#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/*.toml \
           -d 2st.pb_kex_rs \
           -o Output/Single
