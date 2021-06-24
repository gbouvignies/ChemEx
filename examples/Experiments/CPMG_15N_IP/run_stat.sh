#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method_stat.toml \
           -o OutputStat
