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
pathread = 'bin/'
pathreadk = 'bin/klus/'
prct=5e-4
guarda_path = '../pics/'
dpi=400
# ------------------------------------------------------------------------------
# cytwombly_ , groundwater_colorblind_ , contaminant , crazy_mellow , turbo
# tabber , glacier_, gnuplot, batlow, berlin, plasma, prism nipy_spectral
# tokyo, bam, cork, broc roma hawaii lisbon Greys kkrr y-b
colo='abn'
colo_='cuatrocolo' # cuatrocolo qualitcolor

colo_txt = (0,0,0)
fontsize_=12 # 10 12
# ------------------------------------------------------------------------------
#                                 👓👓
# ------------------------------------------------------------------------------
nx = np.fromfile(pathread+'x_size.bin',dtype=int)
nx = nx[0]
nz = np.fromfile(pathread+'z_size.bin',dtype=int)
nz = nz[0]

x = np.fromfile(pathread+'x.bin',dtype=float)
z = np.fromfile(pathread+'z.bin',dtype=float)

psi_ = np.fromfile(pathread+'psi_.bin',dtype=float)
psi_ = np.reshape(psi_,(nx,nz))
psi_ = psi_.T

psi1_ = np.fromfile(pathread+'psi1_.bin',dtype=float)
psi1_ = np.reshape(psi1_,(nx,nz))
psi1_ = psi1_.T

nabmn = np.fromfile(pathread+'pseud_size.bin',dtype=int)
nabmn = nabmn[0]
pseud = np.fromfile(pathread+'pseud.bin',dtype=float)
pseud = np.reshape(pseud,(2,nabmn))
pseud = pseud.T

abmn = np.fromfile(pathread+'abmn.bin',dtype=np.intc)
abmn = np.reshape(abmn,(4,nabmn))
abmn = abmn.T

niabmn = np.fromfile(pathread+'iabmn_size.bin',dtype=int)
niabmn = niabmn[0]
iabmn  = np.fromfile(pathread+'iabmn.bin',dtype=int)
iabmn1 = np.fromfile(pathread+'iabmn1.bin',dtype=int)

nsrc = np.fromfile(pathread+'src_size.bin',dtype=int)
nsrc = nsrc[0]
src  = np.fromfile(pathread+'src.bin',dtype=int)
nrecs= np.fromfile(pathread+'recs_size.bin',dtype=int)
nrecs= nrecs[0]
recs = np.fromfile(pathread+'recs.bin',dtype=int)
recs = np.reshape(recs,(2,nrecs))
nelectrxyz = np.fromfile(pathread+'electrxyz_size.bin',dtype=int)
nelectrxyz= nelectrxyz[0]
electrxyz = np.fromfile(pathread+'electrxyz.bin',dtype=float)
electrxyz = np.reshape(electrxyz,(2,nelectrxyz))
# ------------------------------------------------------------------------------
#                               👓👓 🌂🌂
# ------------------------------------------------------------------------------
nxk = np.fromfile(pathreadk+'x_size.bin',dtype=int)
nxk = nxk[0]
nzk = np.fromfile(pathreadk+'z_size.bin',dtype=int)
nzk = nzk[0]

xk = np.fromfile(pathreadk+'x.bin',dtype=float)
zk = np.fromfile(pathreadk+'z.bin',dtype=float)

psik_ = np.fromfile(pathreadk+'psi_.bin',dtype=float)
psik_ = np.reshape(psik_,(nxk,nzk))
psik_ = psik_.T

nabmnk = np.fromfile(pathreadk+'pseud_size.bin',dtype=int)
nabmnk = nabmnk[0]
pseudk = np.fromfile(pathreadk+'pseud.bin',dtype=float)
pseudk = np.reshape(pseudk,(2,nabmnk))
pseudk = pseudk.T

abmnk = np.fromfile(pathreadk+'abmn.bin',dtype=np.intc)
abmnk = np.reshape(abmnk,(4,nabmnk))
abmnk = abmnk.T

niabmnk = np.fromfile(pathreadk+'iabmn_size.bin',dtype=int)
niabmnk = niabmnk[0]
iabmnk  = np.fromfile(pathreadk+'iabmn.bin',dtype=int)

nsrck = np.fromfile(pathreadk+'src_size.bin',dtype=int)
nsrck = nsrck[0]
srck  = np.fromfile(pathreadk+'src.bin',dtype=int)
nrecsk= np.fromfile(pathreadk+'recs_size.bin',dtype=int)
nrecsk= nrecsk[0]
recsk = np.fromfile(pathreadk+'recs.bin',dtype=int)
recsk = np.reshape(recsk,(2,nrecsk))
nelectrxyzk = np.fromfile(pathreadk+'electrxyz_size.bin',dtype=int)
nelectrxyzk= nelectrxyzk[0]
electrxyzk = np.fromfile(pathreadk+'electrxyz.bin',dtype=float)
electrxyzk = np.reshape(electrxyzk,(2,nelectrxyzk))
# ------------------------------------------------------------------------------
#                                 🚧🎨
# ------------------------------------------------------------------------------
psi_=psi_/prct
psik_=psik_/prct
psi1_=psi1_/prct

midi= 0
mini=-1
maxi= 1

iabmnmin = recs.min()
iabmnmax = recs.max()

iabmnmink = recsk.min()
iabmnmaxk = recsk.max()

extents_xz = [x[0],x[-1],z[-1],z[0]]
extents_xzk = [xk[0],xk[-1],zk[-1],zk[0]]
# ------------------------------------------------------------------------------
nrows = 2
ncols = 3

width_max = x[nx-1]
height_max= z[nz-1]

wi,he = figaspect(nrows*width_max/(ncols*height_max))
# ------------------------------------------------------------------------------
#                                   🎨🎨
#                                   🎨🎨
#                                   🎨🎨
# ------------------------------------------------------------------------------
plt_=fancy_figure(data=np.zeros((2,2)),holdon='close',colorbaron='off').matrix()

fig = plt_.figure()
fig.set_size_inches(wi, he)
ax0 = plt_.subplot2grid((1,3),(0,0),rowspan=1,colspan=1)
ax1 = plt_.subplot2grid((1,3),(0,1),rowspan=1,colspan=1)
ax2 = plt_.subplot2grid((1,3),(0,2),rowspan=1,colspan=1)
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
ax0.annotate("slice in width",xy=(2.6,5.9),fontsize=fontsize_,color=colo_txt,ha='left')
# ------------------------------------------------------------------------------
#  𝐦 pseudo
fancy_figure(
ax_=ax0,
curve=pseud[iabmn1-1,1],
x=pseud[iabmn1-1,0],
symbol='o',
scatter_colo=abmn[iabmn[2]-1,2],
scatter_size=90,
fillstyle='top',
vmin=iabmnmin,
vmax=iabmnmax,
colop=colo_,
colorbaron='off',
holdon='on'
).scatterer()
#  𝐧 pseudo
fancy_figure(
ax_=ax0,
curve=pseud[iabmn1-1,1],
x=pseud[iabmn1-1,0],
symbol='o',
scatter_colo=abmn[iabmn[2]-1,3],
scatter_size=90,
fillstyle='bottom',
vmin=iabmnmin,
vmax=iabmnmax,
colop=colo_,
colorbaron='off',
holdon='on'
).scatterer()
# ------------------------------------------------------------------------------
#                                   🐘 electrxyz 🐘
# ------------------------------------------------------------------------------
# 𝐚
fancy_figure(
ax_=ax0,
curve=electrxyz[1,src[0]-1],
x=electrxyz[0,src[0]-1],
ax_accu="%.0f",
symbol='s',
markersize=6,
colop=((1,1,1)),
holdon='on'
).plotter()
# 𝐛
fancy_figure(
ax_=ax0,
curve=electrxyz[1,src[1]-1],
x=electrxyz[0,src[1]-1],
ax_accu="%.0f",
symbol='s',
markersize=6,
colop=((1,1,1)),
holdon='on'
).plotter()
#  𝐦
fancy_figure(
ax_=ax0,
curve=electrxyz[1,recs[0,0]-1],
x=electrxyz[0,recs[0,0]-1],
symbol='s',
scatter_colo=abmn[iabmn[2]-1,2],
scatter_size=80,
# fillstyle='top',
vmin=iabmnmin,
vmax=iabmnmax,
colo=colo_,
colorbaron='off',
holdon='on'
).scatterer()
#  𝐧
fancy_figure(
ax_=ax0,
curve=electrxyz[1,recs[1,0]-1],
x=electrxyz[0,recs[1,0]-1],
symbol='s',
scatter_colo=abmn[iabmn[2]-1,3],
scatter_size=80,
# fillstyle='bottom',
vmin=iabmnmin,
vmax=iabmnmax,
colop=colo_,
colorbaron='off',
xlabel='Length',
ylabel='Depth',
holdon='on'
).scatterer()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
fancy_figure(
ax_=ax2,
x=xk,
y=zk,
extent=extents_xzk,
data=psik_,
midi=midi,
vmin=mini,
vmax=maxi,
ax_accu="%.0f",
colo=colo,
colorbaron='off',
holdon='on'
).matrix()
# ------------------------------------------------------------------------------
ax2.annotate("slice in width",xy=(0.1,5.9),fontsize=fontsize_,color=colo_txt,ha='left')
# ------------------------------------------------------------------------------
#                                   🐘 electrxyz 🐘
# ------------------------------------------------------------------------------
# 𝐚
fancy_figure(
ax_=ax2,
curve=electrxyzk[1,srck[0]-1],
x=electrxyzk[0,srck[0]-1],
ax_accu="%.0f",
symbol='s',
markersize=6,
colop=((1,1,1)),
holdon='on'
).plotter()
# 𝐛
fancy_figure(
ax_=ax2,
curve=electrxyzk[1,srck[1]-1],
x=electrxyzk[0,srck[1]-1],
ax_accu="%.0f",
symbol='s',
markersize=6,
colop=((1,1,1)),
holdon='on'
).plotter()
#  𝐦
fancy_figure(
ax_=ax2,
curve=electrxyzk[1,recsk[0,:]-1],
x=electrxyzk[0,recsk[0,:]-1],
symbol='s',
scatter_colo=abmnk[iabmnk-1,2],
scatter_size=10*np.ones(iabmnk.shape[0]),
# fillstyle='left',
vmin=iabmnmink,
vmax=iabmnmaxk,
colo=colo_,
colorbaron='off',
holdon='on'
).scatterer()
#  𝐧
fancy_figure(
ax_=ax2,
curve=electrxyzk[1,recsk[1,:]-1],
x=electrxyzk[0,recsk[1,:]-1],
symbol='s',
scatter_colo=abmnk[iabmnk-1,3],
scatter_size=10*np.ones(iabmnk.shape[0]),
# fillstyle='right',
vmin=iabmnmink,
vmax=iabmnmaxk,
colo=colo_,
colorbaron='off',
holdon='on'
).scatterer()
# ------------------------------------------------------------------------------
#  𝐦 pseudolocs
fancy_figure(
ax_=ax2,
curve=pseudk[iabmnk-1,1],
x=pseudk[iabmnk-1,0],
symbol='o',
scatter_colo=abmnk[iabmnk-1,2],
scatter_size=90*np.ones(iabmnk.shape[0]),
fillstyle='top',
vmin=iabmnmink,
vmax=iabmnmaxk,
colo=colo_,
colorbaron='off',
holdon='on'
).scatterer()
#  𝐧 pseudolocs
fancy_figure(
ax_=ax2,
curve=pseudk[iabmnk-1,1],
x=pseudk[iabmnk-1,0],
symbol='o',
scatter_colo=abmnk[iabmnk-1,3],
scatter_size=90*np.ones(iabmnk.shape[0]),
fillstyle='bottom',
vmin=iabmnmink,
vmax=iabmnmaxk,
colo=colo_,
colorbaron='off',
xlabel='Length',
holdon='on'
).scatterer()
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
ax1.annotate("slice in width",xy=(2.6,5.9),fontsize=fontsize_,color=colo_txt,ha='left')
# ------------------------------------------------------------------------------
#                                   🐘 electrxyz 🐘
# ------------------------------------------------------------------------------
# 𝐚
fancy_figure(
ax_=ax1,
curve=electrxyz[1,src[0]-1],
x=electrxyz[0,src[0]-1],
ax_accu="%.0f",
symbol='s',
markersize=6,
colop=((1,1,1)),
holdon='on'
).plotter()
# 𝐛
fancy_figure(
ax_=ax1,
curve=electrxyz[1,src[1]-1],
x=electrxyz[0,src[1]-1],
ax_accu="%.0f",
symbol='s',
markersize=6,
colop=((1,1,1)),
holdon='on'
).plotter()
#  𝐦
fancy_figure(
ax_=ax1,
curve=electrxyz[1,recs[0,:]-1],
x=electrxyz[0,recs[0,:]-1],
symbol='s',
scatter_colo=recs[0,:],
scatter_size=80*np.ones(iabmn.shape[0]),
# fillstyle='top',
vmin=iabmnmin,
vmax=iabmnmax,
colo=colo_,
colorbaron='off',
holdon='on'
).scatterer()
#  𝐧
fancy_figure(
ax_=ax1,
curve=electrxyz[1,recs[1,:]-1],
x=electrxyz[0,recs[1,:]-1],
symbol='s',
scatter_colo=recs[1,:],
scatter_size=80*np.ones(iabmn.shape[0]),
# fillstyle='bottom',
vmin=iabmnmin,
vmax=iabmnmax,
colo=colo_,
colorbaron='off',
holdon='on'
).scatterer()
# ------------------------------------------------------------------------------
#  𝐦 pseudolocs
fancy_figure(
ax_=ax1,
curve=pseud[iabmn-1,1],
x=pseud[iabmn-1,0],
symbol='o',
scatter_colo=abmn[iabmn-1,2],
scatter_size=90*np.ones(iabmn.shape[0]),
fillstyle='top',
vmin=iabmnmin,
vmax=iabmnmax,
colo=colo_,
colorbaron='off',
holdon='on'
).scatterer()
#  𝐧 pseudolocs
fancy_figure(
ax_=ax1,
curve=pseud[iabmn-1,1],
x=pseud[iabmn-1,0],
symbol='o',
scatter_colo=abmn[iabmn-1,3],
scatter_size=90*np.ones(iabmn.shape[0]),
fillstyle='bottom',
vmin=iabmnmin,
vmax=iabmnmax,
colo=colo_,
colorbaron='off',
xlabel='Length',
guarda_path=guarda_path,
guarda=dpi,
fig_name='figo'
).scatterer()
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
im_ = fancy_image(im=im_,nh=nh,nh_r=int(1.9*dpi)).padder_h()
im  = fancy_image(im=im,im_=im_).concat_v()
# im.show()
im.save(guarda_path+"sensis-paper.png","PNG", dpi=dpii)
# ------------------------------------------------------------------------------
print('')
print('')
print('     just finished the 🎨👨🏻‍🎨')
print('')
print('')
# ------------------------------------------------------------------------------
