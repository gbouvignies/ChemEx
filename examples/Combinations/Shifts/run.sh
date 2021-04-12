#!/bin/sh

SELECTION="11 29 85 108 112 113 116 117 134 135 136 139"

chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           --include $SELECTION \
           -o Output
