import scipy.io as sio
import numpy as np
import tensorflow as tf
# ------------------------------------------------------------------------------
print(tf.__version__)
# ------------------------------------------------------------------------------
tf.keras.backend.clear_session()
# ------------------------------------------------------------------------------
ml_name  = 'ml_encodeco'
x_train_ = 'b_mini_'
y_train_ = 'b_mini'
x_test_  = 'c'
# ------------------------------------------------------------------------------
# for alles
path_data = '../../../data/pic2pic/'
path_outi = '../../../data/pic2pic/output/'
# for google-colab
path_data = '/content/'
path_outi = '/content/'
# ------------------------------------------------------------------------------
# 
#                              load data & model
# 
# ------------------------------------------------------------------------------
model = tf.keras.models.load_model(path_outi+ml_name+'.h5')

x_train = sio.loadmat(path_data + x_train_+'.mat')
y_train = sio.loadmat(path_data + y_train_+'.mat')
x_test  = sio.loadmat(path_data + x_test_+'.mat')
x_train = x_train[x_train_]
y_train = y_train[y_train_]
x_test  = x_test[x_test_]
# ------------------------------------------------------------------------------
#
# NOTE:
# 
# [x,y]_train = ( # of elements, size of pictures )
# ------------------------------------------------------------------------------
x_train=np.transpose(x_train,(2,0,1))
y_train=np.transpose(y_train,(2,0,1))
x_test =np.transpose(x_test,(2,0,1))
# ------------------------------------------------------------------------------
n_samples,npix_x,npix_y = x_train.shape
n_tests,_,_ = x_test.shape
# ------------------------------------------------------------------------------
# baby tf cant handle one-channel matrices, 
# so we need to put a fake dimension. 
x_train = np.expand_dims(x_train, axis=-1)
y_train = np.expand_dims(y_train, axis=-1)
x_test  = np.expand_dims(x_test, axis=-1)
# ------------------------------------------------------------------------------
# 
# 
#                   test the model
# 
# 
# ------------------------------------------------------------------------------
#               load all samples & store
# ------------------------------------------------------------------------------
y_reco = np.zeros((n_samples,npix_x,npix_y))
# load and expand dims
for idx in range(0,n_samples):
    x=x_train[idx,:,:,:]
    x=np.expand_dims(x, axis=0)
    # predict
    y=model.predict(x)
    # reduce dims & store
    y_reco[idx,:,:]=y[0,:,:,0]
# ------------------------------------------------------------------------------
#               load all tests & store
# ------------------------------------------------------------------------------
y_test = np.zeros((n_tests,npix_x,npix_y))
# load and expand dims
for idx in range(0,n_tests):
    x=x_test[idx,:,:,:]
    x=np.expand_dims(x, axis=0)
    # predict
    y=model.predict(x)
    # reduce dims & store
    y_test[idx,:,:]=y[0,:,:,0]
# ------------------------------------------------------------------------------
# 
# 
# save everything for easy plotting later
# 
# 
# ------------------------------------------------------------------------------
np.save(path_outi+'y_reco',y_reco)
np.save(path_outi+'y_test',y_test)

np.save(path_outi+'x_train',x_train.squeeze())
np.save(path_outi+'y_train',y_train.squeeze())
np.save(path_outi+'x_test',x_test.squeeze())
# ------------------------------------------------------------------------------
#               load only one sample
# ------------------------------------------------------------------------------
# # load and expand dims
# idx = 70
# x=x_train[idx,:,:,:]
# x=np.expand_dims(x, axis=0)
# # predict
# y=model.predict(x)
# # reduce dims
# y=y[0,:,:,0]
# x=x[0,:,:,0]
# # get true
# y_true=y_train[idx,:,:,0]
# ------------------------------------------------------------------------------
