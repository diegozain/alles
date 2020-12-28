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
# ------------------------------------------------------------------------------
# for alles
path_data = '../../../data/pic2pic/'
path_outi = '../../../data/pic2pic/'
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
# ------------------------------------------------------------------------------
#               load only one sample
# ------------------------------------------------------------------------------
# load and expand dims
idx = 5
x=x_train[idx,:,:,:]
x=np.expand_dims(x, axis=0)
# ------------------------------------------------------------------------------
# 
# 
#                         visualize the journey
# 
# 
# ..............................................................................
# Let's define a new Model that will take an image as input, and will output
# intermediate representations for all layers in the previous model after
# the first.
outputs_     = [layer.output for layer in model.layers[1:]]
model_visual = tf.keras.models.Model(inputs = model.input, outputs = outputs_)
# ..............................................................................
# Let's run our image through our network, thus obtaining all
# intermediate representations for this image.
outputs_ = model_visual.predict(x)
# ..............................................................................
# These are the names of the layers, so can have them as part of our plot
layers_ = [layer.name for layer in model.layers]
# ..............................................................................
# Now let's display our representations
for layer_, output_ in zip(layers_, outputs_):
  # Just do this for the conv / maxpool layers, not the fully-connected layers
  print(output_.shape)
  if len(output_.shape) == 4:
    # The feature map has shape (1, nxy, nxy, nf)
    nf  = output_.shape[-1]
    nxy = output_.shape[1]
    # We will tile our images in this matrix
    display_grid = np.zeros((nxy, nxy * nf))
    # Postprocess the feature to make it visually palatable
    for i_ in range(nf):
      x = output_[0, :, :, i_]
      # x -= x.mean()
      # x /= x.std()
      # x *= 64
      # x += 128
      x = np.clip(x, 0, 255).astype('uint8')
      # We'll tile each filter into this big horizontal grid
      display_grid[:, i_ * nxy : (i_ + 1) * nxy] = x
    # ..........................................................................
    # Display the grid
    scale = 5. / nf
    plt.figure(figsize=(scale * nf, scale))
    plt.title(layer_ , fontsize=20)
    plt.grid(False)
    # ..........................................................................
    ax = plt.gca()
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    # ..........................................................................
    plt.imshow(display_grid, aspect='auto', cmap='plasma')
    # ..........................................................................
    # fig=plt.gcf()
    # fig.savefig(path_outi+layer_+'.png', format='png', dpi=200,bbox_inches='tight',pad_inches=0)
# ..............................................................................