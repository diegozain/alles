import numpy as np
# ------------------------------------------------------------------------------
import sys
sys.path.append('../../../src/python/graphics/fancy_figure/')
from fancy_figure import fancy_figure
sys.path.append('src/')
from bayesian_opti import *
# ------------------------------------------------------------------------------
#
#
# this is a script that shows the method with an example.
# 
# in most cases, we cannot evaluate the objective function 
# inside the bayesian optimization scheme.
# 
# if that is the case, you should use the function:
# 
#     x = optimize_baye(x,y,bounds)
# 
# y are the values of the objective function
# x are the parameters that were used to compute each value of y
# bounds are the bounds for each parameter
#
# x is of size samples by parameters
# y is of size samples by one
# bounds is of size parameters by two ( two because of [min, max])
# 
#
#
# ------------------------------------------------------------------------------
np.random.seed(seed=0)
size=[6,4]
# ------------------------------------------------------------------------------
# 
# 
#               set up true objective function 
# 
# 
# 
# ------------------------------------------------------------------------------
# complicated example
# ------------------------------------------------------------------------------
def obj_fnc(x):
    y = ((x[:,0]-2)**2+(x[:,1]-2)**2 - 0.5) * ((x[:,0]+2)**2+(x[:,1]+1.5)**2-0.05) * ((x[:,0]-1.5)**2 + (x[:,1]+2)**2 + 0.2)
    y = y.reshape(-1,1)
    return y
def obj_fnc_(xx,yy):
    y = ((xx-2)**2+(yy-2)**2 - 0.5) * ((xx+2)**2+(yy+1.5)**2-0.05)*((xx-1.5)**2 + (yy+2)**2 + 0.2) 
    return y

midi=200
# ------------------------------------------------------------------------------
# very simple example
# ------------------------------------------------------------------------------
# def obj_fnc(x):
#     y = ((x[:,0]-2)**2+(x[:,1]-2)**2)
#     return y.reshape(-1,1)
# def obj_fnc_(xx,yy):
#     return ((xx-2)**2+(yy-2)**2)
# 
# midi=5
# ------------------------------------------------------------------------------
# 
# 
#               set up bayesian opti 
# 
# 
# 
# ------------------------------------------------------------------------------
n_init=5
x_samples = np.random.uniform(-5, 5, size=(n_init, 2))
y_samples = obj_fnc(x_samples)
# ------------------------------------------------------------------------------
# iterations
n_iter = 5
bounds = np.array([[-5, 5],[-5, 5]])
# ------------------------------------------------------------------------------
x_samples,y_samples = optimize_baye_(x_samples,y_samples,obj_fnc,bounds,n_iter)
# ------------------------------------------------------------------------------
print('found arg min = ',x_samples[-1,:])
print('found minimum = ',-y_samples[-1])
# ------------------------------------------------------------------------------
# 
# 
# 
#                           for plotting only
# 
# 
# 
# ------------------------------------------------------------------------------
# for plotting true objective function
# ------------------------------------------------------------------------------
x_param = np.arange(-5, 5, 0.1)
y_param = np.arange(-5, 5, 0.1)
# ------------------------------------------------------------------------------
xx,yy = np.meshgrid(x_param, y_param, sparse=True)
y_true = obj_fnc_(xx,yy)

print('true minimum  = ',y_true.min())
# ------------------------------------------------------------------------------
# for plotting approximated (last) objective function
# ------------------------------------------------------------------------------
# Gaussian process
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, RBF
ker_ = ConstantKernel(1.0) * RBF(length_scale=1)
gau_pro_ = GaussianProcessRegressor(kernel=ker_,alpha=1e-0)
gau_pro_.fit(x_samples, y_samples)

x_approx_ = np.arange(-5, 5, 0.1)
y_approx_ = np.arange(-5, 5, 0.1)
nx=x_approx_.size
ny=y_approx_.size
y_approx = np.zeros([nx,ny])
for ix in range(nx):
    for iy in range(ny):
        tmp_ = np.array([x_approx_[ix],y_approx_[iy]])
        tmp_ = tmp_.reshape(1, -1)
        tmp_ = exp_imp(tmp_, x_samples, y_samples, gau_pro_)
        y_approx[iy,ix] = tmp_[0]
# ------------------------------------------------------------------------------
# 
# 
# 
#                   actual plotting
# 
# 
# 
# 
# ------------------------------------------------------------------------------
plt_=fancy_figure(
data=np.zeros((2,2)),holdon='close',colorbaron='off').matrix()
# ------------------------------------------------------------------------------
fig,ax=plt_.subplots(1,2)#,figsize=size)
# ------------------------------------------------------------------------------
extents_ = [-5,5,5,-5]

fancy_figure(
ax_=ax[0],
x=2.06,
curve=2.14,
symbol='*',
colop='b',
markersize=10,
holdon='on'
).plotter()

fancy_figure(
ax_=ax[0],
x=x_samples[0:n_init,0],
curve=x_samples[0:n_init,1],
symbol='.',
colop='y',
markersize=13,
holdon='on'
).plotter()

fancy_figure(
ax_=ax[0],
x=x_samples[n_init:-1,0],
curve=x_samples[n_init:-1,1],
symbol='.',
colop='g',
markersize=8,
holdon='on'
).plotter()

fancy_figure(
ax_=ax[0],
x=x_samples[-1,0],
curve=x_samples[-1,1],
symbol='*',
colop='r',
markersize=8,
holdon='on'
).plotter()

fancy_figure(
ax_=ax[0],
extent=extents_,
data=y_true,
x=x_param,
y=y_param,
midi=midi,
colo='Greys',
colorbaron='off',
xlabel='Parameter \#1',
ylabel='Parameter \#2',
title='True',
holdon='on'
).matrix()
# ------------------------------------------------------------------------------
fancy_figure(
ax_=ax[1],
x=2.06,
curve=2.14,
symbol='*',
colop='b',
markersize=10,
holdon='on'
).plotter()

fancy_figure(
ax_=ax[1],
x=x_samples[0:n_init,0],
curve=x_samples[0:n_init,1],
symbol='.',
colop='y',
markersize=8,
holdon='on'
).plotter()

fancy_figure(
ax_=ax[1],
x=x_samples[n_init:-1,0],
curve=x_samples[n_init:-1,1],
symbol='.',
colop='g',
markersize=8,
holdon='on'
).plotter()

fancy_figure(
ax_=ax[1],
x=x_samples[-1,0],
curve=x_samples[-1,1],
symbol='*',
colop='r',
markersize=8,
holdon='on'
).plotter()

fancy_figure(
figsize=size,
ax_=ax[1],
extent=extents_,
data=-y_approx,
x=x_approx_,
y=y_approx_,
colo='Greys',
colorbaron='off',
xlabel='Parameter \#1',
y_ticklabels='off',
title='Approximate',
fig_name='bayes-opti-ex',
guarda_path='../pics/',
guarda=200
).matrix()
# ------------------------------------------------------------------------------
