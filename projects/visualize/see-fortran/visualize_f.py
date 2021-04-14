import sys
sys.path.append('../../../src/python/graphics/fancy_figure/')
from fancy_figure import fancy_figure
# ------------------------------------------------------------------------------
size=[4*2,4]
guarda_path = '../pics/'
dpi=120

nrows = 3
ncols = 2
# ------------------------------------------------------------------------------
import numpy as np
# ------------------------------------------------------------------------------
# # load when saved as characters
# t = np.loadtxt('t.dat')
# one_d= np.loadtxt('one_d.dat')
# two_d= np.loadtxt('two_d.dat')
# 
# two_d=two_d.reshape((ncols,nrows))
# two_d=two_d.transpose()
# ------------------------------------------------------------------------------
# load when saved as binary
t = np.fromfile('t.dat',dtype=float)
one_d = np.fromfile('one_d.dat',dtype=float)
two_d = np.fromfile('two_d.dat',dtype='int32')

two_d=two_d.reshape((ncols,nrows))
two_d=two_d.transpose()
# ------------------------------------------------------------------------------
print('shape of 1D data ',one_d.shape)
print('shape of 2D data ',two_d.shape)
# ------------------------------------------------------------------------------
# 
#                    plotting routine begins 
# 
# ------------------------------------------------------------------------------
# 2D matrix
fancy_figure(
aspect='auto',
figsize=[5,5],
data=two_d,
x_ticklabels='off',
y_ticklabels='off',
title='2D data',
ylabel='\# of rows',
xlabel='\# of columns',
guarda_path=guarda_path,
guarda=dpi,
fig_name='fortran-2D'
).matrix()
# ------------------------------------------------------------------------------
fancy_figure(
figsize=size,
curve=one_d,
x=t,
title='1D data',
xlabel='Time (time units)',
ylabel='Amplitude',
symbol='-',
colo_accu='%2i',
colop=((0.0471,0.4824,0.8627)),
margin=(0.05,0.1),
guarda_path=guarda_path,
guarda=dpi,
fig_name='fortran-1D'
).plotter()
# ------------------------------------------------------------------------------
