#!/bin/sh

SELECTION="3 4 12 13 14 30 32 65 68 70 72 74 75 78 80 91 100 101 103 104 105 108 110 111 113 115 116 117 118 119 121 122 137 138 139 141 142 146 149"

chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -d 2st_rs \
           --include $SELECTION \
           -o Output
