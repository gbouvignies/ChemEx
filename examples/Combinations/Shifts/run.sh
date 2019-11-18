#!/bin/sh

chemex fit -e Experiments/cp*.toml Experiments/sqmq*.toml \
           -p Parameters/*.toml \
           --include 11 29 85 108 112 113 116 117 134 135 136 139 \
           -o Output
