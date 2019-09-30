#!/bin/bash

rm -rf Output/Single

for i in 3 4 7 8 9 11 15 25 26 27 29 30 31 32 35 47 48 62 65 67 68 69 70 71 72 73 74 77 78 80 81 82 84 90 91 94 95 96 97 98 99 104 108 110 111 113 114 115 116 117 118 119 120 123 124 125 126 133 135 137 139 140 141 142 145 146 148 150 151 152 153 154 155 158 162 164
do
chemex fit \
    -e Experiments/*.toml \
    -p Parameters/*.toml \
    -m Methods/method.toml \
    -r $i \
    -o Output/Single/${i}H
done

