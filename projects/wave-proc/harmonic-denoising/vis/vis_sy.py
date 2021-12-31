import sys
sys.path.append('../../../../src/python/graphics/fancy_figure/')
from fancy_figure import fancy_figure
# ------------------------------------------------------------------------------
size=[4*2,4]
guarda_path = ''
dpi=120
# ------------------------------------------------------------------------------
import numpy as np
# ------------------------------------------------------------------------------
# load when saved as binary
uo_obs = np.fromfile('../bin/uo_obs.bin',dtype=float)
uo_rec= np.fromfile('../bin/uo_rec.bin',dtype=float)
uo_sig= np.fromfile('../bin/uo_sig.bin',dtype=float)
t = np.fromfile('../bin/t.bin',dtype=float)
# ------------------------------------------------------------------------------
#
#                    plotting routine begins
#
# ------------------------------------------------------------------------------
fancy_figure(
figsize=size,
curve=uo_obs,
x=t,
symbol='-',
# colop=((0.5,0,0)),
margin=(0.05,0.1),
holdon = 'on'
).plotter()

fancy_figure(
figsize=size,
curve=uo_sig,
x=t,
symbol='-',
# colop=((0.5,0,0)),
margin=(0.05,0.1),
holdon = 'on'
).plotter()

fancy_figure(
figsize=size,
curve=uo_rec,
x=t,
title='Harmonic denoising',
xlabel='Time (s)',
ylabel='Amplitude',
symbol='--',
# colop=((0.5,0,0)),
margin=(0.05,0.1),
# guarda_path=guarda_path,
# guarda=dpi,
# fig_name='fortran-1D'
).plotter()
# ------------------------------------------------------------------------------
