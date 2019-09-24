#!/bin/sh

for i in 14 16 18 35 36 43 51 54 57 60
do
  chemex fit -e Experiments/*.toml \
             -p Parameters/global.toml \
                Parameters/cs_a.toml \
                Parameters/dw_ab.toml \
             -m Methods/method.toml \
             -o Output/${i}H-C \
             -r $i
done
