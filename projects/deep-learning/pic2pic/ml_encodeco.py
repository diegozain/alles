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
#                       optimization parameters
# 
# ------------------------------------------------------------------------------
step_=0.1
n_epoch = 24 #256 24
epochs_ = 250
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
n_samples,npix_x,npix_y = x_train.shape
# ------------------------------------------------------------------------------
# baby tf cant handle one-channel matrices, 
# so we need to put a fake dimension. 
x_train = np.expand_dims(x_train, axis=-1)
y_train = np.expand_dims(y_train, axis=-1)

x_input = (npix_y,npix_x,1)
# ------------------------------------------------------------------------------
# 
#                     the model (learning machine)
# 
# ------------------------------------------------------------------------------
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
model  = generator(x_input)
# ------------------------------------------------------------------------------
# gives a nice overview of what the model is like
model.summary()
# ------------------------------------------------------------------------------
# optimizer = tf.keras.optimizers.SGD(lr=step_)
# ------------------------------------------------------------------------------
# for Adadelta
# learning_rate is the gradient step-size
# rho is the decay rate (?)
# epsilon is a regularization constant
# -------------------------------
optimizer = tf.keras.optimizers.Adadelta(
    learning_rate=step_, rho=0.95, epsilon=1e-07
    )
# ------------------------------------------------------------------------------
model.compile(optimizer=optimizer,
              metrics  =['accuracy'],
              loss     ='MeanSquaredLogarithmicError'#'MeanSquaredError' 
              )
# ------------------------------------------------------------------------------
history = model.fit(x=x_train, y=y_train,
                    steps_per_epoch=n_samples,
                    epochs = epochs_,
                    batch_size=1,
                    verbose=0)
# ------------------------------------------------------------------------------
loss = history.history['loss']
loss = np.asarray(loss)
accu = history.history['accuracy']
accu = np.asarray(accu)
# ------------------------------------------------------------------------------
np.save(path_outi+'loss',loss)
np.save(path_outi+'accu',accu)
model.save(path_outi+'ml_encodeco.h5')
# ------------------------------------------------------------------------------