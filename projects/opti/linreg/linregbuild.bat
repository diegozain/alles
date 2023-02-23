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
REM ECHO "       ***************  compiling  ****************"
REM ----------------------------------------------------------------------------
del *.mod *.obj *.exe
REM ----------------------------------------------------------------------------
ifort /Qmkl /Qopenmp /c /traceback /heap-arrays "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\lapack.f90" linreg.f90 linreg_.f90 parafit.f90 linefit.f90

ifort linreg.obj lapack.obj
ifort linreg_.obj lapack.obj
ifort parafit.obj lapack.obj
ifort linefit.obj lapack.obj

ifort /Qmkl /Qopenmp /c /traceback /heap-arrays ..\..\..\src\fortran\qrfits.f90 qrfits_ie.f90

ifort qrfits_ie.obj lapack.obj qrfits.obj
REM ----------------------------------------------------------------------------
del *.mod *.obj
REM ----------------------------------------------------------------------------
REM ECHO "       ******** author: diego domenzain  **********"
REM ----------------------------------------------------------------------------
REM .\linreg.exe
REM .\linreg_.exe
REM .\parafit.exe
REM ----------------------------------------------------------------------------
