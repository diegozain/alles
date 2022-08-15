@echo off
setlocal enableDelayedExpansion
REM ----------------------------------------------------------------------------
REM
REM                i think windows is 💩💩💩💩
REM
REM   diego 2022
REM ----------------------------------------------------------------------------
set myvar=somevar
set mytxt=..\..\
echo %PATH%
REM ----------------------------------------------------------------------------
REM for /l %%x in (1, 1, 2) do (
REM    echo %%x
REM    set file=%myvar%%%x.txt
REM    type nul >> !file!
REM    REM >> appends
REM    REM >  overwrites
REM    echo !mytxt!>> !file!
REM    REM echo !file!> !file!
REM )
REM ----------------------------------------------------------------------------
set clus=1 2 4 5
for %%x in (%clus%) do (
  echo doing clus %%x
  set file=%myvar%%%x.txt
  type nul >> !file!
  REM >> appends
  REM >  overwrites
  echo !mytxt!>> !file!
  REM echo !file!> !file!
)
REM ----------------------------------------------------------------------------
exit /B
REM ----------------------------------------------------------------------------
