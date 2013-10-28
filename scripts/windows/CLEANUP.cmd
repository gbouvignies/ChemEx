@echo off

if exist chemex_pkg. (
  echo. && echo --- Deleting build --- && echo.
  rd /Q /S build
  echo. && echo --- Deleting install files --- && echo.
  del /Q %HOMEPATH%\AppData\Roaming\Python\Python27\site-packages\chemex*
  del /Q %HOMEPATH%\AppData\Roaming\Python\Scripts\chemex*
  echo. && echo --- Deleting chemex_dumps from examples --- && echo.
  for /D /R %%X in (chemex_dump*) do rd /S /Q %%X
  echo. && echo --- Deleting .pyc files --- && echo.
  del /S /Q *.pyc > nul
  echo. && echo --- Attempting ChemEx ... --- && echo.
  chemex_fit.py
  timeout 5
) else (
  echo. && echo --- You are in the wrong directory! Please return to src! && echo.
)
