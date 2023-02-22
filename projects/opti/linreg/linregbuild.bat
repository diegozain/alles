@ECHO OFF
REM ----------------------------------------------------------------------------
REM
REM                ifort in 💩
REM
REM   diego 2021
REM ----------------------------------------------------------------------------
REM C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\mkl_pardiso.f90
REM cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
REM ----------------------------------------------------------------------------
REM del *.mod *.obj *.exe
REM ----------------------------------------------------------------------------
ifort /Qmkl /Qopenmp /c /traceback /heap-arrays "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\lapack.f90" linreg.f90 linreg_.f90 parafit.f90

ifort linreg.obj lapack.obj
ifort linreg_.obj lapack.obj
ifort parafit.obj lapack.obj
REM ----------------------------------------------------------------------------
REM .\linreg.exe
REM .\linreg_.exe
REM .\parafit.exe
REM ----------------------------------------------------------------------------
