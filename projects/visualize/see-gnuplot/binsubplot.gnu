#!/usr/bin/gnuplot
# -------------------------------------------------------------------------
#                                matlab
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
# 
# in2mm = 25.4 # mm (fixed)
# pt2mm = 0.3528 # mm (fixed)
# 
# mm2px = dpi/in2mm
# ptscale = pt2mm*mm2px
# round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
# wpx = round(width * mm2px)
# hpx = round(height * mm2px)
# 
# set terminal pngcairo size wpx,hpx fontscale ptscale linewidth ptscale pointscale ptscale
# -------------------------------------------------------------------------
# # # wxt
# # set terminal wxt size 700,524 enhanced font 'Verdana,20' persist
# png
set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output '../pics/binsubplot.png'
# -------------------------------------------------------------------------
unset key
# -------------------------------------------------------------------------
# 
#                                macros
# 
# -------------------------------------------------------------------------
set macros
POS = "at graph 0.92,0.9 font ',8'"
# borders & style
STYLE="set style line 11 lc rgb '#808080' lt 1"
BORDER="set border 3 front ls 11"
TICKS="set tics nomirror out scale 0.75"
# axis image & tight
IMAGE="set size ratio -1"
TIGHT="set autoscale fix"
# -------------------------------------------------------------------------
# 
#                                🎨🎨🎨
# 
# -------------------------------------------------------------------------
set multiplot layout 1,2 title "😀𝙿Ω" font ",14"
# -------------------------------------------------------------------------
unset key
# decor
@STYLE; @BORDER; @TICKS; @IMAGE; @TIGHT
# labels
set xlabel 'x (µm)'
set ylabel 'z (µm)'
set title '🥸' 
set title font ",15"
unset colorbox
# not cluttered axis
stats "srcgnuplot/x.bin" binary format="%float" u 1 nooutput
xmin = ceil(STATS_min)
xmax = floor(STATS_max)
set xtics xmin , ((xmax-xmin)/4.0) , xmax
stats "srcgnuplot/z.bin" binary format="%float" u 1 nooutput
zmin = ceil(STATS_min)
zmax = floor(STATS_max)
set ytics zmin , ((zmax-zmin)/4.0) , zmax
# subplot label
set label front 'a' at xmin*1.01,zmax font ',20' textcolor rgb '#ffffff'
# plot
plot 'srcgnuplot/sigxz.bin' binary matrix with image
# -------------------------------------------------------------------------
unset key
# decor
@STYLE; @BORDER; @TICKS; @IMAGE; @TIGHT
# labels
unset colorbox
unset ylabel
# not cluttered axis
stats "srcgnuplot/x.bin" binary format="%float" u 1 nooutput
xmin = ceil(STATS_min)
xmax = floor(STATS_max)
set xtics xmin , ((xmax-xmin)/4.0) , xmax
stats "srcgnuplot/z.bin" binary format="%float" u 1 nooutput
zmin = ceil(STATS_min)
zmax = floor(STATS_max)
set ytics zmin , ((zmax-zmin)/4.0) , zmax
# plot
plot 'srcgnuplot/sigxz.bin' binary matrix with image
# -------------------------------------------------------------------------
unset multiplot
# -------------------------------------------------------------------------
