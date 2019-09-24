#!/bin/sh

rm -rf output_single

for i in 1 2 5 6 7 8 11 14 15 17 19 21 26 32 36 38 39 43 44 45 46 52 53 61 62 63 64 65 67 68 69 71 72 73 74 77 78 80 83 85 87 88
do
  chemex fit \
    -e Experiments/*.toml \
    -p Parameters/*.toml \
    -m Methods/method_single.toml \
    -o Output/Single/${i}N-H \
    -r $i
done
