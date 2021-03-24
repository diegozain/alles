#!/bin/bash

gfortran gen_data.f90
./a.out
python visualize_f.py

rm *.o *.mod *.out *.dat fort.*
