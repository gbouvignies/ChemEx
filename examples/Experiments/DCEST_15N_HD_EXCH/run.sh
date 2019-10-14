#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/*.toml \
           -m Methods/*.toml \
           -d 2st.hd_exch \
           -o Output
