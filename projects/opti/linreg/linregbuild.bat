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
ifort /Qmkl /c "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\lapack.f90" linreg.f90
ifort /Qmkl linreg.obj lapack.obj

ifort /Qmkl /Qopenmp /c linreg_.f90
ifort /Qmkl /Qopenmp linreg_.obj lapack.obj

ifort /Qmkl /c parafit.f90
ifort /Qmkl parafit.obj lapack.obj
REM ----------------------------------------------------------------------------
REM .\linreg.exe
REM .\linreg_.exe
REM .\parafit.exe
REM ----------------------------------------------------------------------------
