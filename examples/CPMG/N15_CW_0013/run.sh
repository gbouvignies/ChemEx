#!/bin/bash

chemex fit \
    -e Experiments/*.cfg \
    -p Parameters/params_n15.cfg \
    -m Methods/method_global.cfg \
    -o Output
