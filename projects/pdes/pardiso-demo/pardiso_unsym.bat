@ECHO OFF
REM ----------------------------------------------------------------------------
REM
REM                ifort
REM
REM   diego 2021
REM ----------------------------------------------------------------------------
del *.mod *.obj *.exe
ifort /Qmkl pardiso_unsym.f90
.\pardiso_unsym.exe
del *.mod *.obj
REM ----------------------------------------------------------------------------