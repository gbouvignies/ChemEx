#!/bin/sh

chemex fit -e Experiments/19Hz_500ms_*cfg \
           -p Parameters/par.cfg \
           -m Methods/met.cfg \
           -d 2st.pb_kex \
           -o Output
