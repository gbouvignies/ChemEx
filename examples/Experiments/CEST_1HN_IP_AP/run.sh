#!/bin/sh

for i in 3 4 12 13 14 30 32 65 68 70 72 74 75 78 80 91 100 101 103 104 105 108 110 111 113 115 116 117 118 119 121 122 137 138 139 141 142 146 149
do
  chemex fit -e Experiments/*.toml \
             -p Parameters/global.toml \
                Parameters/cs_a.toml \
                Parameters/dw_ab.toml \
             -m Methods/method.toml \
             -o Output/${i}HN-N \
             -r ${i}HN-N
done
