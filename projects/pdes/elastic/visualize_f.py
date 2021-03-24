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
sz_time= np.loadtxt('sz_time.dat')
rx= np.loadtxt('rx.dat')
vz_rx= np.loadtxt('vz_rx.dat')
vx_rx= np.loadtxt('vx_rx.dat')
# ------------------------------------------------------------------------------
rx = rx[:,0]
nt = t.size
nrx= rx.size

vx_rx_real = vx_rx[:,0]
vx_rx_imag = vx_rx[:,1]
vx_rx_real=vx_rx_real.reshape((nrx,nt))
vx_rx_imag=vx_rx_imag.reshape((nrx,nt))
vx_rx_real=vx_rx_real.transpose()
vx_rx_imag=vx_rx_imag.transpose()

vz_rx_real = vz_rx[:,0]
vz_rx_imag = vz_rx[:,1]
vz_rx_real=vz_rx_real.reshape((nrx,nt))
vz_rx_imag=vz_rx_imag.reshape((nrx,nt))
vz_rx_real=vz_rx_real.transpose()
vz_rx_imag=vz_rx_imag.transpose()
# ------------------------------------------------------------------------------
# 
#                    plotting routine begins 
# 
# ------------------------------------------------------------------------------
fancy_figure(
figsize=size,
curve=sz_time[:,0],
x=t,
title=r'Source on $v_z$',
xlabel='Time (s)',
ylabel='Amplitude',
symbol='-',
colop=((0.3,0.3,0.3)),
margin=(0.05,0.1),
guarda_path=guarda_path,
guarda=dpi,
fig_name='elastic-lamb-source'
).plotter()
# ------------------------------------------------------------------------------
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
data=vx_rx_real,
x=rx,
y=t,
midi=0,
# colo_accu='%1.e',
colorbaron='off',
title=r'Receivers $v_x$',
ylabel='Time (s)',
xlabel='Length from source (m)',
holdon='on'
).matrix()

fancy_figure(
ax_=ax[1],
aspect='auto',
extent=extent_,
figsize=[8,5],
data=vz_rx_real,
x=rx,
y=t,
midi=0,
# colo_accu='%1.e',
colorbaron='off',
title=r'Receivers $v_z$',
# ylabel='Time (s)',
xlabel='Length from source (m)',
guarda_path=guarda_path,
guarda=dpi,
fig_name='elastic-lamb-data'
).matrix()
# ------------------------------------------------------------------------------
# 
# 
# 
# ------------------------------------------------------------------------------
plt_=fancy_figure(data=np.zeros((2,2)),holdon='close',colorbaron='off').matrix()
fig,ax=plt_.subplots(1,2,figsize=size)
# ------------------------------------------------------------------------------
vx_rx_real=vx_rx_real/vx_rx_real.max()
drx=rx[1]-rx[0]
for irx in range(0,nrx-1):
    fancy_figure(
    ax_=ax[0],
    figsize=[8,5],
    curve = 0.5*drx*vx_rx_real[:,irx]+rx[irx],
    # curve=vx_rx_real[:,irx] + irx,
    x=t,
    symbol='-',
    colo_accu='%2i',
    colop=((0.3,0.3,0.3)),
    holdon='on'
    ).plotter()
    
fancy_figure(
ax_=ax[0],
figsize=[8,5],
curve = 0.5*drx*vx_rx_real[:,nrx-1]+rx[nrx-1],
# curve=vx_rx_real[:,nrx-1] + nrx-1,
x=t,
title=r'Receivers $v_x$',
xlabel='Time (s)',
ylabel='Length from source (m)',
symbol='-',
colo_accu='%2i',
colop=((0.3,0.3,0.3)),
margin=(0.05,0.1),
holdon='on'
).plotter()
# ------------------------------------------------------------------------------
vz_rx_real=vz_rx_real/vz_rx_real.max()
drx=rx[1]-rx[0]
for irx in range(0,nrx-1):
    fancy_figure(
    ax_=ax[1],
    figsize=[8,5],
    curve = 0.5*drx*vz_rx_real[:,irx]+rx[irx],
    # curve=vz_rx_real[:,irx] + irx,
    x=t,
    symbol='-',
    colo_accu='%2i',
    colop=((0.3,0.3,0.3)),
    holdon='on'
    ).plotter()
    
fancy_figure(
ax_=ax[1],
figsize=[8,5],
curve = 0.5*drx*vz_rx_real[:,nrx-1]+rx[nrx-1],
# curve=vz_rx_real[:,nrx-1] + nrx-1,
x=t,
title=r'Receivers $v_z$',
xlabel='Time (s)',
# ylabel='Length from source (m)',
symbol='-',
colo_accu='%2i',
colop=((0.3,0.3,0.3)),
margin=(0.05,0.1),
guarda_path=guarda_path,
guarda=dpi,
fig_name='elastic-lamb-data-'
).plotter()
# ------------------------------------------------------------------------------
