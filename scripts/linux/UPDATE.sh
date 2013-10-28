#!/bin/bash
hg pull
hg update
read -t5 -n1 -r -p "Press any key. Continuing in five seconds..." key
./INSTALL
