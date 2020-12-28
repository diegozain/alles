import numpy as np
# ------------------------------------------------------------------------------
import sys
sys.path.append('../../../src/python/graphics/fancy_figure/')
from fancy_figure import fancy_figure
# ------------------------------------------------------------------------------
#
#
# this script assumes you already have all outputs.
# 
# it only plots things beautifully. omg how cute.
# 
# 
# 
#
#
# ------------------------------------------------------------------------------
path_outi = '../../../data/pic2pic/output/'
# ------------------------------------------------------------------------------
accu = np.load(path_outi+'accu'+'.npy')
loss = np.load(path_outi+'loss'+'.npy')

# input outside the training data
x_test = np.load(path_outi+'x_test'+'.npy')
# input of the training data
x_train= np.load(path_outi+'x_train'+'.npy')

# true output of the training data
y_train= np.load(path_outi+'y_train'+'.npy')
# recovered output of the training data
y_reco = np.load(path_outi+'y_reco'+'.npy')
# recovered output outside of the training data
y_test = np.load(path_outi+'y_test'+'.npy')
# ------------------------------------------------------------------------------
maxi = max([x_test.max(),x_train.max(),y_train.max(),y_test.max(),y_reco.max()])
mini = min([x_test.min(),x_train.min(),y_train.min(),y_test.min(),y_reco.min()])

n_samples,npix_x,npix_y = x_train.shape
n_tests,_,_ = x_test.shape
# ------------------------------------------------------------------------------
plt_=fancy_figure(
data=np.zeros((2,2)),holdon='close',colorbaron='off').matrix()
# ------------------------------------------------------------------------------
fig,ax=plt_.subplots(4,6)#,figsize=size)
# ------------------------------------------------------------------------------
for i_ in range(n_samples-1):
	ic = int(i_%6)
	ir = int((i_-ic)/6)
	
	fancy_figure(
	ax_=ax[ir,ic],
	data=y_train[i_,:,:],
	vmin=mini,
	vmax=maxi,
	x_ticklabels='off',
	y_ticklabels='off',
	# colo='Greys',
	colorbaron='off',
	holdon='on'
	).matrix()
	
i_=i_+1
ic = int(i_%6)
ir = int((i_-ic)/6)

fancy_figure(
ax_=ax[ir,ic],
data=y_train[i_,:,:],
vmin=mini,
vmax=maxi,
x_ticklabels='off',
y_ticklabels='off',
# colo='Greys',
colorbaron='off',
super_title='True (trained)',
fig_name='pic2pic_true',
guarda_path='../pics/',
guarda=200
).matrix()
# ------------------------------------------------------------------------------
fig,ax=plt_.subplots(4,6)#,figsize=size)
# ------------------------------------------------------------------------------
for i_ in range(n_samples-1):
	ic = int(i_%6)
	ir = int((i_-ic)/6)
	
	fancy_figure(
	ax_=ax[ir,ic],
	data=y_reco[i_,:,:],
	vmin=mini,
	vmax=maxi,
	x_ticklabels='off',
	y_ticklabels='off',
	# colo='Greys',
	colorbaron='off',
	holdon='on'
	).matrix()
	
i_=i_+1
ic = int(i_%6)
ir = int((i_-ic)/6)

fancy_figure(
ax_=ax[ir,ic],
data=y_reco[i_,:,:],
vmin=mini,
vmax=maxi,
x_ticklabels='off',
y_ticklabels='off',
# colo='Greys',
colorbaron='off',
super_title='Recovered (trained)',
fig_name='pic2pic_reco',
guarda_path='../pics/',
guarda=200
).matrix()
# ------------------------------------------------------------------------------
fig,ax=plt_.subplots(4,6)#,figsize=size)
# ------------------------------------------------------------------------------
for i_ in range(n_samples-1):
	ic = int(i_%6)
	ir = int((i_-ic)/6)
	
	fancy_figure(
	ax_=ax[ir,ic],
	data=x_train[i_,:,:],
	vmin=mini,
	vmax=maxi,
	x_ticklabels='off',
	y_ticklabels='off',
	# colo='Greys',
	colorbaron='off',
	holdon='on'
	).matrix()
	
i_=i_+1
ic = int(i_%6)
ir = int((i_-ic)/6)

fancy_figure(
ax_=ax[ir,ic],
data=x_train[i_,:,:],
vmin=mini,
vmax=maxi,
x_ticklabels='off',
y_ticklabels='off',
# colo='Greys',
colorbaron='off',
super_title='Observed (trained)',
fig_name='pic2pic_obse',
guarda_path='../pics/',
guarda=200
).matrix()
# ------------------------------------------------------------------------------
fig,ax=plt_.subplots(1,5)#,figsize=size)
# ------------------------------------------------------------------------------
for i_ in range(n_tests-1):
	
	fancy_figure(
	ax_=ax[i_],
	data=x_test[i_,:,:],
	vmin=mini,
	vmax=maxi,
	x_ticklabels='off',
	y_ticklabels='off',
	# colo='Greys',
	colorbaron='off',
	holdon='on'
	).matrix()
	
i_=i_+1

fancy_figure(
ax_=ax[i_],
data=x_test[i_,:,:],
vmin=mini,
vmax=maxi,
x_ticklabels='off',
y_ticklabels='off',
# colo='Greys',
colorbaron='off',
super_title='Observed (not trained)',
fig_name='pic2pic_test_x',
guarda_path='../pics/',
guarda=200
).matrix()
# ------------------------------------------------------------------------------
fig,ax=plt_.subplots(1,5)#,figsize=size)
# ------------------------------------------------------------------------------
for i_ in range(n_tests-1):
	
	fancy_figure(
	ax_=ax[i_],
	data=y_test[i_,:,:],
	vmin=mini,
	vmax=maxi,
	x_ticklabels='off',
	y_ticklabels='off',
	# colo='Greys',
	colorbaron='off',
	holdon='on'
	).matrix()
	
i_=i_+1

fancy_figure(
ax_=ax[i_],
data=y_test[i_,:,:],
vmin=mini,
vmax=maxi,
x_ticklabels='off',
y_ticklabels='off',
# colo='Greys',
colorbaron='off',
super_title='Recovered (not trained)',
fig_name='pic2pic_test_y',
guarda_path='../pics/',
guarda=200
).matrix()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
