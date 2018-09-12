#!/bin/sh

chemex fit -e Experiments/*.cfg \
           -p Parameters/params.cfg \
           -d 2st.pb_kex \
           -o Output
