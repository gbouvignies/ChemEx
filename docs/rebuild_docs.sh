#!/bin/sh

rm -rf build/
rm -rf source/modules/*/_autosummary/
make html
make latexpdf
make epub
#make man
#make info
