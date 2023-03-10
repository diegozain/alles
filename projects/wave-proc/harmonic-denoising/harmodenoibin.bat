@ECHO OFF
REM ----------------------------------------------------------------------------
REM
REM                ifort
REM
REM   diego 2021
REM ----------------------------------------------------------------------------
REM C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\mkl_pardiso.f90
REM cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'

REM del *.mod *.obj *.exe
del *.mod *.obj

REM /traceback /check:bounds

ifort /Qmkl /Qopenmp /c /traceback /heap-arrays ..\..\..\src\fortran\calculus.f90 ..\..\..\src\fortran\readfiles.f90 ..\..\..\src\fortran\harmodenoiser.f90 harmodenoibin.f90

ifort harmodenoibin.obj calculus.obj readfiles.obj harmodenoiser.obj

REM .\harmodenoibin.exe

REM del *.mod *.obj *.exe
del *.mod *.obj
REM ----------------------------------------------------------------------------
