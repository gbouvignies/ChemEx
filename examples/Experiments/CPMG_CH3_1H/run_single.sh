#!/bin/bash

for i in 1HD1-CD1 2HD1-CD1 3HD1-CD1 4HD1-CD1 6HD1-CD1 7HD1-CD1 9HD1-CD1 10HD1-CD1 14HD1-CD1 15HD1-CD1
do
  chemex fit \
    -e Experiments/*.toml \
    -p Parameters/global.toml \
       Parameters/cs_a.toml \
       Parameters/dw_ab.toml \
    -m Methods/method_global.toml \
    -o Output/Single/"$i" \
    +r "$i"
done

