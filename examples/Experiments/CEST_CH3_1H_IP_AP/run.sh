#!/bin/sh

for i in L7HD1-CD1 L7HD2-CD2 I29HD1-CD1 V71HG1-CG1 I100HD1-CD1 M102HE-CE V103HG1-CG1 M106HE-CE L118HD2-CD2 L121HD1-CD1 L121HD2-CD2
do
  chemex fit -e Experiments/*.toml \
             -p Parameters/global.toml \
                Parameters/cs_a.toml \
                Parameters/dw_ab.toml \
                Parameters/j_a.toml \
             -m Methods/method.toml \
             -o Output/$i \
             -r $i
done
