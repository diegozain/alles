import sys
sys.path.append('../../../src/python/graphics/fancy_figure/')
sys.path.append('../../../src/python/graphics/fancy_image/')
from fancy_figure import fancy_figure
from fancy_image import fancy_image
import numpy as np
import numpy.matlib
from math import nan
from matplotlib.figure import figaspect
# ------------------------------------------------------------------------------
path_read = 'bin/'
prct=5e-4
guarda_path = '../pics/'
dpi=350
# ------------------------------------------------------------------------------
# cytwombly_ , groundwater_colorblind_ , contaminant , crazy_mellow , turbo
# tabber , glacier_, gnuplot, batlow, berlin, plasma, tokyo, bam
colo='cork'

colo_txt = (0,0,0)
fontsize_=12 # 10 12
# ------------------------------------------------------------------------------
#                                 👓👓
# ------------------------------------------------------------------------------
nx = np.fromfile(path_read+'x_size.bin',dtype=int)
nx = nx[0]
nz = np.fromfile(path_read+'z_size.bin',dtype=int)
nz = nz[0]

x = np.fromfile(path_read+'x.bin',dtype=float)
z = np.fromfile(path_read+'z.bin',dtype=float)

psi_ = np.fromfile(path_read+'psi_.bin',dtype=float)
psi_ = np.reshape(psi_,(nx,nz))
psi_ = psi_.T

psi1_ = np.fromfile(path_read+'psi1_.bin',dtype=float)
psi1_ = np.reshape(psi1_,(nx,nz))
psi1_ = psi1_.T

nabmn = np.fromfile(path_read+'pseud_size.bin',dtype=int)
nabmn = nabmn[0]
pseud = np.fromfile(path_read+'pseud.bin',dtype=float)
pseud = np.reshape(pseud,(2,nabmn))
pseud = pseud.T

niabmn = np.fromfile(path_read+'iabmn_size.bin',dtype=int)
niabmn = niabmn[0]
iabmn = np.fromfile(path_read+'iabmn.bin',dtype=int)
iabmn1 = np.fromfile(path_read+'iabmn1.bin',dtype=int)
# ------------------------------------------------------------------------------
#                                 🚧🎨
# ------------------------------------------------------------------------------
psi_=psi_/prct
psi1_=psi1_/prct

midi= 0
mini=-1
maxi= 1

extents_xz = [x[0],x[-1],z[-1],z[0]]
# ------------------------------------------------------------------------------
nrows = 2
ncols = 3

width_max = x[nx-1]
height_max= z[nz-1]

wi,he = figaspect(nrows*width_max/(ncols*height_max))
# ------------------------------------------------------------------------------
#                                 🎨🎨
# ------------------------------------------------------------------------------
plt_=fancy_figure(data=np.zeros((2,2)),holdon='close',colorbaron='off').matrix()

fig = plt_.figure()
fig.set_size_inches(wi, he)
ax0 = plt_.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
ax1 = plt_.subplot2grid((1,2),(0,1),rowspan=1,colspan=1)
# ------------------------------------------------------------------------------
# fig, ax = plt_.subplots()
# fig.suptitle(title_mega,fontsize=20)
# ------------------------------------------------------------------------------
fancy_figure(
ax_=ax0,
x=x,
y=z,
extent=extents_xz,
data=psi1_,
midi=midi,
vmin=mini,
vmax=maxi,
ax_accu="%.0f",
colo=colo,
colorbaron='off',
holdon='on'
).matrix()
# ------------------------------------------------------------------------------
ax0.annotate("slice in width",xy=(3.2,5.9),fontsize=fontsize_,color=colo_txt,ha='left')
# ------------------------------------------------------------------------------
#  pseudolocs
fancy_figure(
ax_=ax0,
curve=pseud[iabmn1-1,1],
x=pseud[iabmn1-1,0],
ax_accu="%.0f",
symbol='.',
markersize=20,
colop=((0.4940,0.1840,0.5560)),
holdon='on'
).plotter()
#
fancy_figure(
ax_=ax0,
curve=pseud[iabmn1-1,1],
x=pseud[iabmn1-1,0],
ax_accu="%.0f",
symbol='.',
markersize=10,
colop=((0.9290,0.6940,0.1250)),
xlabel='Length',
ylabel='Depth',
holdon='on'
).plotter()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
fancy_figure(
ax_=ax1,
x=x,
y=z,
extent=extents_xz,
data=psi_,
midi=midi,
vmin=mini,
vmax=maxi,
ax_accu="%.0f",
colo=colo,
colorbaron='off',
holdon='on'
).matrix()
# ------------------------------------------------------------------------------
ax1.annotate("slice in width",xy=(3.2,5.9),fontsize=fontsize_,color=colo_txt,ha='left')
# ------------------------------------------------------------------------------
#  pseudolocs
fancy_figure(
ax_=ax1,
curve=pseud[iabmn-1,1],
x=pseud[iabmn-1,0],
ax_accu="%.0f",
symbol='.',
markersize=20,
colop=((0.4940,0.1840,0.5560)),
holdon='on'
).plotter()
#
fancy_figure(
ax_=ax1,
curve=pseud[iabmn-1,1],
x=pseud[iabmn-1,0],
ax_accu="%.0f",
symbol='.',
markersize=10,
colop=((0.9290,0.6940,0.1250)),
xlabel='Length',
# ylabel='Depth',
guarda_path=guarda_path,
guarda=dpi,
fig_name='figo'
).plotter()
# ------------------------------------------------------------------------------
fig = plt_.gcf()
plt_.close()
size = fig.get_size_inches()
fancy_figure(
midi=midi,
vmin=mini,
vmax=maxi,
colo=colo,
# colo_accu="%.0f",
colo_title='Normalized sensitivity',
figsize=size,
holdon='on',
guarda_path=guarda_path,
guarda=dpi).colorbar_alone()
# ------------------------------------------------------------------------------
# join figures 👺
# ------------------------------------------------------------------------------
im =guarda_path+"figo.png"
im_=guarda_path+"fig.png"
im = fancy_image(im=im).openim()
dpii=im.info['dpi']
nh,_=im.size
im_ = fancy_image(im=im_,nh=nh,nh_r=int(1.5*dpi)).padder_h()
im  = fancy_image(im=im,im_=im_).concat_v()
# im.show()
im.save(guarda_path+"sensis-paper.png","PNG", dpi=dpii)
# ------------------------------------------------------------------------------
