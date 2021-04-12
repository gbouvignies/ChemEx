#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -d 2st_hd \
           -o Output
