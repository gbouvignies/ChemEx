#!/bin/sh

chemex simulate -e Experiments/*.toml \
                -p Parameters/parameters.toml \
                -o OutputSim
