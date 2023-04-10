#!/usr/bin/gnuplot
# -------------------------------------------------------------------------
reset
# -------------------------------------------------------------------------
# # wxt
# set terminal wxt size 700,524 enhanced font 'Verdana,20' persist
# png
set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
set output '../pics/binplot.png'
# -------------------------------------------------------------------------
unset key

# Axes
set tics nomirror out scale 0.75
# Grid
set style line 12 lc rgb'#808080' lt 0 lw 1
set grid back ls 12

# line color
set style line 1 lt 1 lw 2 pt 3 linecolor rgb '#e94d1d'
set size square

set xlabel 'x (µm)'
set ylabel 'y (µm)'
# -------------------------------------------------------------------------
# # this plots an array with index coords in the x axis
# plot "uo_sig.bin" binary format="%double" using 0:1 w l ls 1
# -------------------------------------------------------------------------
stats "srcgnuplot/uost.bin" binary format="%double" u 1 nooutput
x_lim = ceil(STATS_max)
set xrange[-0.2: x_lim+0.2]
# set xtics x_lim/4.0
stats "srcgnuplot/uost.bin" binary format="%double" u 2 nooutput
y_lim = ceil(STATS_max)
set yrange[-0.2: y_lim+0.2]
# set ytics y_lim/4.0

# this plots xy axis from uost, but uost has to be saved as (2,nt) !!
# gnuplot follows c convention, not fortran
plot "srcgnuplot/uost.bin" binary format="%double" using 1:2 w l ls 1
# -------------------------------------------------------------------------



