#!/bin/sh

chemex fit -e Experiments/*.cfg \
           -p Parameters/*.cfg \
           -m Methods/*.cfg \
           -d 2st.pb_kex \
           -o Output
