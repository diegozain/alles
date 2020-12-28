import scipy.io as sio
import numpy as np
import tensorflow as tf
# ------------------------------------------------------------------------------
print(tf.__version__)
# ------------------------------------------------------------------------------
tf.keras.backend.clear_session()
tf.random.set_seed(51)
np.random.seed(51)
# ------------------------------------------------------------------------------
x_train_ = 'b_mini_'
y_train_ = 'b_mini'
# ------------------------------------------------------------------------------
path_data = '../../../data/pic2pic/'
path_outi = '../../../data/pic2pic/'
path_data = '/content/'
path_outi = '/content/'
# ------------------------------------------------------------------------------
# 
#                              load data
# 
# ------------------------------------------------------------------------------
x_train = sio.loadmat(path_data + x_train_+'.mat')
y_train = sio.loadmat(path_data + y_train_+'.mat')
x_train = x_train[x_train_]
y_train = y_train[y_train_]
# ------------------------------------------------------------------------------
#
# NOTE:
# 
# [x,y]_train = ( # of elements, size of pictures )
# ------------------------------------------------------------------------------
x_train=np.transpose(x_train,(2,0,1))
y_train=np.transpose(y_train,(2,0,1))
# ------------------------------------------------------------------------------
n_samples,npix_y,npix_x = x_train.shape
# ------------------------------------------------------------------------------
# baby tf cant handle one-channel matrices, 
# so we need to put a fake dimension. 
x_train = np.expand_dims(x_train, axis=-1)
y_train = np.expand_dims(y_train, axis=-1)
# ------------------------------------------------------------------------------
# 
#                       optimization parameters
# 
# ------------------------------------------------------------------------------
nfx=3
nfy=3

n_epoch = 24
n_iter = 250

x_input = (npix_y,npix_x,1)
y_output= (npix_y,npix_x,1)
n_batch = 3
# ------------------------------------------------------------------------------
# 
# 
# 
#                     the model (learning machine)
# 
# 
# 
# ------------------------------------------------------------------------------
# 
# 						discriminator model
# 
# ------------------------------------------------------------------------------
# define the standalone discriminator model
def discriminator(x_input,y_output):
	# weight initialization
	weight_init = tf.keras.initializers.RandomNormal(stddev=0.02)
	# --------------------------------------------------------------------------
	# inputs
	# --------------------------------------------------------------------------
	# x input
	x_input = tf.keras.layers.Input(shape=x_input)
	# y output
	y_output= tf.keras.layers.Input(shape=y_output)
	# concat input and output
	x_y_together = tf.keras.layers.Concatenate()([x_input, y_output])
	# --------------------------------------------------------------------------
	# convolutions
	# --------------------------------------------------------------------------
	classi = tf.keras.layers.Conv2D(2, (4,4), strides=(1,1),
	 padding='same', kernel_initializer=weight_init)(x_y_together)
	classi = tf.keras.layers.LeakyReLU(alpha=0.2)(classi)
	#
	classi = tf.keras.layers.Conv2D(4, (4,4), strides=(1,1),
	 padding='same', kernel_initializer=weight_init)(classi)
	classi = tf.keras.layers.BatchNormalization()(classi)
	classi = tf.keras.layers.LeakyReLU(alpha=0.2)(classi)
	#
	classi = tf.keras.layers.Conv2D(8, (4,4), strides=(1,1),
	 padding='same', kernel_initializer=weight_init)(classi)
	classi = tf.keras.layers.BatchNormalization()(classi)
	classi = tf.keras.layers.LeakyReLU(alpha=0.2)(classi)
	#
	classi = tf.keras.layers.Conv2D(16, (4,4), strides=(1,1),
	 padding='same', kernel_initializer=weight_init)(classi)
	classi = tf.keras.layers.BatchNormalization()(classi)
	classi = tf.keras.layers.LeakyReLU(alpha=0.2)(classi)
	#
	classi = tf.keras.layers.Conv2D(1, (4,4), padding='same', kernel_initializer=weight_init,activation='sigmoid')(classi)
	# --------------------------------------------------------------------------
	# ml model
	# --------------------------------------------------------------------------
	# define model
	d_model = tf.keras.models.Model([y_output,x_input], classi)
	# define optimizer
	optimizer_ = tf.keras.optimizers.Adadelta(learning_rate=0.01, rho=0.95, epsilon=1e-07)
	# compile model
	d_model.compile(loss='binary_crossentropy', optimizer=optimizer_, metrics=['accuracy'])
	return d_model
# ------------------------------------------------------------------------------
# 
# 						generator model
# 
# ------------------------------------------------------------------------------
# define the standalone generator model
def generator(x_input):
	# weight initialization
	weight_init = tf.keras.initializers.RandomNormal(stddev=0.02)
	# --------------------------------------------------------------------------
	# input
	# --------------------------------------------------------------------------
	# x_input input
	x_input = tf.keras.layers.Input(shape=x_input)
	# --------------------------------------------------------------------------
	# convolutions (encoder)
	# --------------------------------------------------------------------------
	enco_1 = tf.keras.layers.Conv2D(2, (4,4), strides=(1,1),
   padding='same', kernel_initializer=weight_init)(x_input)
	enco_1 = tf.keras.layers.LeakyReLU(alpha=0.2)(enco_1)
	#
	enco_2 = tf.keras.layers.Conv2D(4, (4,4), strides=(1,1),
   padding='same', kernel_initializer=weight_init)(enco_1)
	enco_2 = tf.keras.layers.BatchNormalization()(enco_2)
	enco_2 = tf.keras.layers.LeakyReLU(alpha=0.2)(enco_2)
	#
	enco_3 = tf.keras.layers.Conv2D(8, (4,4), strides=(1,1),
   padding='same', kernel_initializer=weight_init)(enco_2)
	enco_3 = tf.keras.layers.BatchNormalization()(enco_3)
	enco_3 = tf.keras.layers.LeakyReLU(alpha=0.2)(enco_3)
	#
	enco_4 = tf.keras.layers.Conv2D(16, (4,4), strides=(1,1),
   padding='same', kernel_initializer=weight_init)(enco_3)
	enco_4 = tf.keras.layers.BatchNormalization()(enco_4)
	enco_4 = tf.keras.layers.LeakyReLU(alpha=0.2)(enco_4)
	# --------------------------------------------------------------------------
	# convolutions (bottleneck)
	# --------------------------------------------------------------------------
	x_enco = tf.keras.layers.Conv2D(16, (4,4), strides=(1,1),
   padding='same', kernel_initializer=weight_init,activation='relu')(enco_4)
	x_enco = tf.keras.layers.Activation('relu')(x_enco)
	# --------------------------------------------------------------------------
	# convolutions (decoder)
	# --------------------------------------------------------------------------
	deco_1 = tf.keras.layers.Conv2DTranspose(16, (4,4), strides=(1,1),
   padding='same', kernel_initializer=weight_init)(x_enco)
	deco_1 = tf.keras.layers.BatchNormalization()(deco_1)
	deco_1 = tf.keras.layers.Concatenate()([deco_1, enco_4])
	deco_1 = tf.keras.layers.Activation('relu')(deco_1)
	#
	deco_2 = tf.keras.layers.Conv2DTranspose(8, (4,4), strides=(1,1),
   padding='same', kernel_initializer=weight_init)(deco_1)
	deco_2 = tf.keras.layers.BatchNormalization()(deco_2)
	deco_2 = tf.keras.layers.Concatenate()([deco_2, enco_3])
	deco_2 = tf.keras.layers.Activation('relu')(deco_2)
	#
	deco_3 = tf.keras.layers.Conv2DTranspose(4, (4,4), strides=(1,1),
   padding='same', kernel_initializer=weight_init)(deco_2)
	deco_3 = tf.keras.layers.BatchNormalization()(deco_3)
	deco_3 = tf.keras.layers.Concatenate()([deco_3, enco_2])
	deco_3 = tf.keras.layers.Activation('relu')(deco_3)
	#
	deco_4 = tf.keras.layers.Conv2DTranspose(2, (4,4), strides=(1,1),
   padding='same', kernel_initializer=weight_init)(deco_3)
	deco_4 = tf.keras.layers.BatchNormalization()(deco_4)
	deco_4 = tf.keras.layers.Concatenate()([deco_4, enco_1])
	deco_4 = tf.keras.layers.Activation('relu')(deco_4)
	# --------------------------------------------------------------------------
	# output y
	# --------------------------------------------------------------------------
	y_output = tf.keras.layers.Conv2DTranspose(1, (4,4), strides=(1,1),
   padding='same', kernel_initializer=weight_init,activation='tanh')(deco_4)
	# --------------------------------------------------------------------------
	# ml model
	# --------------------------------------------------------------------------
	# define model
	g_model = tf.keras.models.Model(x_input, y_output)
	return g_model
# ------------------------------------------------------------------------------
# 
# 						generator + discriminator model
# 
# ------------------------------------------------------------------------------
# define the combined generator and discriminator model, for updating the generator
def gan_(g_model, d_model):
	# make weights in the discriminator not trainable
	d_model.trainable = False
	
	x_input = g_model.input
	y_output= g_model.output
	classi  = d_model.output
	
	# define gan model as taking noise and label and outputting a classification
	gan_model = tf.keras.models.Model(x_input, [classi, y_output])
	# optimizer
	optimizer_= tf.keras.optimizers.Adadelta(learning_rate=0.01, rho=0.95, epsilon=1e-07)
	# compile model
	gan_model.compile(loss=['binary_crossentropy', 'mae'], optimizer=optimizer_,loss_weights=[1,100])
	return gan_model
# ------------------------------------------------------------------------------
# 
# 						generate samples
# 
# ------------------------------------------------------------------------------
# select true examples
def reals_(x_train,y_train,n_batch):
	# choose random instances
	i_label = np.random.randint(0, y_train.shape[0], n_batch,n_patch)
	# select y_train and x_train
	y_train, x_train = y_train[i_label,:,:,:], x_train[i_label,:,:,:]
	# generate class x_train
	classi = np.ones((n_batch, n_patch, n_patch, 1))
	return [x_train, y_train], classi
# ------------------------------------------------------------------------------
# use the generator to generate n fake examples, with class labels
def fakes_(g_model,x_train,n_batch,n_patch):
	# predict outputs
	y_train = g_model.predict(x_train)
	# create class labels
	classi = np.zeros((n_batch, n_patch, n_patch, 1))
	return y_train, classi
# ------------------------------------------------------------------------------
# 
# 						train
# 
# ------------------------------------------------------------------------------
# train the generator and discriminator
def train(g_model,d_model,gan_model,x_train,y_train,n_iter,n_batch):
	n_patch = d_model.output_shape[0]
	g_loss = np.zeros(n_iter)
	# iterations
	for i_ in range(n_iter):
		# get randomly selected 'real' samples
		[x_train_, y_train_], classi = reals_(x_train,y_train,n_batch,n_patch)
		# update discriminator model weights
		d_loss1 = d_model.train_on_batch([x_train_,y_train_], classi)
		# generate 'fake' examples
		y_train_, classi = fakes_(g_model,x_train_,n_batch,n_patch)
		# update discriminator model weights
		d_loss2 = d_model.train_on_batch([x_train_,y_train_], classi)
		# update the generator via the discriminator's error
		g_loss_,_,_ = gan_model.train_on_batch(x_train_,[classi,y_train_])
		# summarize performance
		print('>%d, d1[%.3f] d2[%.3f] g[%.3f]' % (i_, d_loss1, d_loss2, g_loss_))
		# save loss
		g_loss[i_]= g_loss_
	# save the models
	g_model.save(path_outi+'g_model.h5')
	d_model.save(path_outi+'d_model.h5')
	gan_model.save(path_outi+'gan_model.h5')
	# save the loss
	np.save(path_outi+'g_loss',g_loss)
# ------------------------------------------------------------------------------
d_model  = discriminator(x_input, y_output)
g_model  = generator(x_input)
gan_model= gan_(g_model, d_model)
# ------------------------------------------------------------------------------
# print nice things
d_model.summary()
g_model.summary()
gan_model.summary()
# ------------------------------------------------------------------------------
# train model
train(g_model,d_model,gan_model,x_train,y_train,n_iter,n_batch)
# ------------------------------------------------------------------------------

idx=20
x=x_train[idx,:,:,:]
x=np.expand_dims(x, axis=0)

y=g_model.predict(x)
y=y[0,:,:,0]
x=x[0,:,:,0]

y_true=y_train[idx,:,:,0]

import matplotlib.pyplot as plt

plt.figure()
plt.imshow(x, vmin=0, vmax=1, cmap='jet',)
plt.colorbar()
plt.title('Input' , fontsize=20)

plt.figure()
plt.imshow(y, cmap='jet', vmin=0, vmax=1)
plt.colorbar()
plt.title('Output' , fontsize=20)

plt.figure()
plt.imshow(y_true, vmin=0, vmax=1, cmap='jet',)
plt.colorbar()
plt.title('True' , fontsize=20)

plt.figure()
plt.plot(np.log10(loss))
plt.title('Loss' , fontsize=20)
print('\n last iteration\n',loss[-1])