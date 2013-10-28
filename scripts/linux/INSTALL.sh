#!/bin/bash
echo "Installing ChemEx ... (see install.log)"
python setup.py install --user > install.log 2> errors.log
if [ `wc errors.log | awk '{print $1}'` != 0 ] 
then
  echo "Install had errors ..."
  echo -e "\033[31m" `cat errors.log` "\033[0m"
else
  echo
  echo "RUNNING unittests"
  echo
  sleep 1
  chemex_test.py 2> test.log
  egrep "=====|FAILED|failures|errors|Traceback" -A 50 test.log 
fi
