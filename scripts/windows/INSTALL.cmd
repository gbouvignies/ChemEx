@echo off

python setup.py install --user
echo. && echo ---- NOW  RUNNING  unittests ---- && echo.
timeout 3
chemex_test.py
pause
