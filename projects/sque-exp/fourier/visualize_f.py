import sys
sys.path.append('../../../src/python/graphics/fancy_figure/')
from fancy_figure import fancy_figure
# ------------------------------------------------------------------------------
size=[4*2,4]
guarda_path = '../pics/'
dpi=120
# ------------------------------------------------------------------------------
import numpy as np
f = np.loadtxt('frequency.dat')
t = np.loadtxt('time.dat')
g_= np.loadtxt('data_f.dat')
g = np.loadtxt('data_t.dat')
# ------------------------------------------------------------------------------
# 
#                    plotting routine begins 
# 
# ------------------------------------------------------------------------------
plt_=fancy_figure(data=np.zeros((2,2)),holdon='close',colorbaron='off').matrix()
# ------------------------------------------------------------------------------
fig,ax=plt_.subplots(1,2,figsize=size)
# ------------------------------------------------------------------------------
fancy_figure(
ax_=ax[0],
figsize=size,
# curve=g_[:,0],
curve=np.fft.fftshift(g_[:,0]),
x=f,
title=r'$\mathbf{R}$',
xlabel='Frequency (Hz)',
ylabel='Amplitude',
symbol='-',
colo_accu='%2i',
colop=((0.0471,0.4824,0.8627)),
margin=(0.1,0.1),
holdon='on'
).plotter()

fancy_figure(
ax_=ax[1],
figsize=size,
curve=np.fft.fftshift(g_[:,1]),
x=f,
title=r'$i\mathbf{R}$',
xlabel='Frequency (Hz)',
symbol='-',
colo_accu='%2i',
colop=((0.8627,0.1961,0.1255)),
margin=(0.1,0.1),
guarda_path=guarda_path,
guarda=dpi,
fig_name='data-f'
).plotter()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
fancy_figure(
figsize=size,
curve=g[:,0],
x=t,
title='Time domain',
xlabel='Time (s)',
ylabel='Amplitude',
symbol='-',
colo_accu='%2i',
colop=((0.3,0.3,0.3)),
margin=(0.01,0.1),
guarda_path=guarda_path,
guarda=dpi,
fig_name='data-t'
).plotter()
# ------------------------------------------------------------------------------
