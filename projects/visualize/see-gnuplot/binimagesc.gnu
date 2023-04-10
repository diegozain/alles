#!/usr/bin/gnuplot
# -------------------------------------------------------------------------
# matlab output
# 
# x=(0:0.1:3*pi);
# z=(0:0.1:5*pi);
# [xx,zz]=meshgrid(x,z);
# sig=sin(xx).*cos(zz);
# 
# figure;
# imagesc(sig);
# axis image
# 
# save_bin('sig',sig.','double');
#
# nx=numel(x);
# nz=numel(z);
# sigxz=zeros(nx+1,nz+1);
# sigxz(1,1)=nx;
# sigxz(1,2:nz+1)=z;
# sigxz(2:nx+1,1)=x.';
# sigxz(2:nx+1,2:nz+1)=sig.';
# save_bin('sigxz',sigxz,'float');
# save_bin('x',x,'float');
# save_bin('z',z,'float');
# -------------------------------------------------------------------------
reset
# -------------------------------------------------------------------------
# dpi = 300 ## dpi (variable)
# width = 700 ## mm (variable)
# height = 524 ## mm (variable)

# in2mm = 25.4 # mm (fixed)
# pt2mm = 0.3528 # mm (fixed)

# mm2px = dpi/in2mm
# ptscale = pt2mm*mm2px
# round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
# wpx = round(width * mm2px)
# hpx = round(height * mm2px)

# set terminal pngcairo size wpx,hpx fontscale ptscale linewidth ptscale pointscale ptscale
# -------------------------------------------------------------------------
# # # wxt
# # set terminal wxt size 700,524 enhanced font 'Verdana,20' persist
# png
set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
# -------------------------------------------------------------------------
set output '../pics/binimagesc.png'
# -------------------------------------------------------------------------

unset key

# border
set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror out scale 0.75

set size ratio -1
set autoscale fix


# colorbar
# disable colorbar tics
# set cbtics scale 0

set xlabel 'x (µm)'
set ylabel 'z (µm)'

set title '🥸' 
set title font ",30"

set cblabel "Units (sec²)"
# -------------------------------------------------------------------------
# # classic imagesc
# # index numbers for x&y axis
# plot "srcgnuplot/sig.bin" binary array=(95,158) format="%double" with image
# -------------------------------------------------------------------------
# x&y are stored in sigxz.
# sigxz has to be saved as "float" idk why
# so all are floats here.

stats "srcgnuplot/x.bin" binary format="%float" u 1 nooutput
xmin = ceil(STATS_min)
xmax = floor(STATS_max)
set xtics xmin , ((xmax-xmin)/4.0) , xmax
stats "srcgnuplot/z.bin" binary format="%float" u 1 nooutput
zmin = ceil(STATS_min)
zmax = floor(STATS_max)
set ytics zmin , ((zmax-zmin)/4.0) , zmax

plot 'srcgnuplot/sigxz.bin' binary matrix with image
# -------------------------------------------------------------------------
# multiplot
# set multiplot layout 2,2 title "😀ℵ" font ",14"
# unset key
# plot 'srcgnuplot/sigxz.bin' binary matrix with image
# unset key
# plot 'srcgnuplot/sigxz.bin' binary matrix with image
# unset key
# plot 'srcgnuplot/sigxz.bin' binary matrix with image
# unset key
# plot 'srcgnuplot/sigxz.bin' binary matrix with image
# -------------------------------------------------------------------------


