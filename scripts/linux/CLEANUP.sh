#!/bin/bash
if [ -e chemex ]
then
 echo "Removing \'build\'"
 rm -fr build
 echo "Removing install files"
 rm -fr ~/.local/lib/python2.7/site-packages/chemex*
 rm -fr ~/.local/bin/chemex_*
 echo "Removing chemex_dumps from examples"
 rm -fr `find . -name "chemex_dump*"`
 echo "Removing .pyc files"
 rm -f `find chemex/ -name "*.pyc"`
 echo "Attempting ChemEx ..."
# chemex_fit.py
else 
 echo "You are in the wrong directory! Please go the the main directory"
 echo "(containing the README file) and run the script from there."
fi


