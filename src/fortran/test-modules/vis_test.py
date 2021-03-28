import sys
sys.path.append('../../../src/python/graphics/fancy_figure/')
from fancy_figure import fancy_figure
# ------------------------------------------------------------------------------
size=[4*2,4]
guarda_path = ''
dpi=120
# ------------------------------------------------------------------------------
import numpy as np
t = np.loadtxt('t.dat')
y = np.loadtxt('y.dat')
y_diff= np.loadtxt('y_diff.dat')
y_inte= np.loadtxt('y_inte.dat')
# ------------------------------------------------------------------------------
nt = t.size
# ------------------------------------------------------------------------------
# # here you transpose/reshape vector to matrix
# data=data.reshape((ncols,nrows))
# data=data.transpose()
# ------------------------------------------------------------------------------
# 
#                    plotting routine begins 
# 
# ------------------------------------------------------------------------------
#                             differentiate 
fancy_figure(
figsize=size,
curve=y,
x=t,
symbol='-',
colop=((0.3,0.5,0.3)),
holdon='on'
).plotter()

fancy_figure(
figsize=size,
curve=y_diff,
x=t,
title=r'Differentiate',
xlabel='Time',
ylabel='Amplitude',
symbol='-',
colop=((0.5,0.3,0.3)),
margin=(0.05,0.1),
guarda_path=guarda_path,
guarda=dpi,
fig_name='differentiate'
).plotter()
# ------------------------------------------------------------------------------
#                             integrate 
fancy_figure(
figsize=size,
curve=y,
x=t,
symbol='-',
colop=((0.3,0.5,0.3)),
holdon='on'
).plotter()

fancy_figure(
figsize=size,
curve=y_inte,
x=t,
title=r'Integrate',
xlabel='Time',
ylabel='Amplitude',
symbol='-',
colop=((0.5,0.3,0.3)),
margin=(0.05,0.1),
guarda_path=guarda_path,
guarda=dpi,
fig_name='integrate'
).plotter()
# ------------------------------------------------------------------------------
'''
extent_ = np.array([rx[0],rx[-1],t[-1],t[0]],dtype=float)
# ------------------------------------------------------------------------------
plt_=fancy_figure(data=np.zeros((2,2)),holdon='close',colorbaron='off').matrix()
fig,ax=plt_.subplots(1,2,figsize=size)
# ------------------------------------------------------------------------------
fancy_figure(
ax_=ax[0],
aspect='auto',
extent=extent_,
figsize=[8,5],
data=data,
x=x,
y=y,
midi=0,
# colo_accu='%1.e',
colorbaron='off',
title=r'Data',
ylabel='Columns',
xlabel='Rows',
holdon='on'
).matrix()

fancy_figure(
ax_=ax[1],
aspect='auto',
extent=extent_,
figsize=[8,5],
data=data,
x=x,
y=y,
midi=0,
# colo_accu='%1.e',
colorbaron='off',
title=r'Data',
# ylabel='Time (s)',
xlabel='Columns',
guarda_path=guarda_path,
guarda=dpi,
fig_name='matrix'
).matrix()
'''
# ------------------------------------------------------------------------------