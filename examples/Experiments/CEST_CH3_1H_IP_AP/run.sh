#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/*.toml \
           -m Methods/*.toml \
           -d 2st.pb_kex_rs \
           -o Output
