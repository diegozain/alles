import numpy as np
from matplotlib.colors import ListedColormap
# diego domenzain
# 12.2024
# ---------------------------------------------------------------------------------
# import sys
# sys.path.append("../src/")
# from coloritobar import coloritobar
# colorito = coloritobar()
# 
# coloritomap = colorito.build(-1,0,1,255)
# ---------------------------------------------------------------------------------
class coloritobar:
  def __init__(self):
    pass

  def build(self, minval, midval, maxval, nrgb, cmap="214"):
    # ------------------------------------------------------------------------------
    #
    #                                    🔴🔶🔷
    #                                 build colormap
    #                                       σ
    # 
    # h="B4FBB8"
    # a=tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
    # ------------------------------------------------------------------------------
    if cmap=="214":
      # 🟣🔴 214
      coloritosr = np.asarray([9/255,160/255,0.9,113/255,176/255])
      coloritosg = np.asarray([9/255,136/255,0.9,133/255,90/255])
      coloritosb = np.asarray([11/255,17/255,0.9,145/255,58/255])
    elif cmap=="77":
      # 🟣🔴 77
      coloritosr = np.asarray([203/255,112/255,40/255])
      coloritosg = np.asarray([209/255,35/255,150/255])
      coloritosb = np.asarray([39/255,163/255,224/255])
    elif cmap=="124":
      # 🟣🔴 124
      coloritosr = np.asarray([139/255,154/255,232/255,86/255])
      coloritosg = np.asarray([237/255,234/255,25/255,55/255])
      coloritosb = np.asarray([230/255,100/255,63/255,40/255])
    elif cmap=="128":
      # 🟣🔴 128
      coloritosr = np.asarray([37/255,87/255,244/255,43/255])
      coloritosg = np.asarray([141/255,163/255,74/255,37/255])
      coloritosb = np.asarray([206/255,83/255,22/255,35/255])
    elif cmap=="152":
      # 🟣🔴 152
      coloritosr = np.asarray([20/255,88/255,232/255,232/255,75/255,247/255,42/255,142/255])
      coloritosg = np.asarray([32/255,155/255,134/255,206/255,228/255,203/255,119/255,79/255])
      coloritosb = np.asarray([145/255,95/255,232/255,134/255,242/255,44/255,37/255,28/255])
    elif cmap=="bi":
      coloritosr = np.asarray([116/255,45/255])
      coloritosg = np.asarray([221/255,90/255])
      coloritosb = np.asarray([95/255,226/255])
    else:
      raise ValueError("that color map is not available.")
    # ------------------------------------------------------------------------------
    xx = np.linspace(0,1,coloritosr.shape[0])
    # ------------------------------------------------------------------------------
    r_ = np.interp(np.linspace(0,1,nrgb),xx,coloritosr)
    g_ = np.interp(np.linspace(0,1,nrgb),xx,coloritosg)
    b_ = np.interp(np.linspace(0,1,nrgb),xx,coloritosb)
    # ------------------------------------------------------------------------------
    h = np.asarray([minval , midval , maxval])
    xx = np.asarray([0 , 0.5 , 1])
    yy = np.linspace(0,1,nrgb) 

    f = np.interp(yy,xx,h)

    gg = np.linspace(minval,maxval,nrgb)
    gfofx = np.interp(gg,f,yy)

    r = np.interp(gfofx,yy,r_)
    g = np.interp(gfofx,yy,g_)
    b = np.interp(gfofx,yy,b_)
    # ------------------------------------------------------------------------------
    #                           map colormap
    # ------------------------------------------------------------------------------
    coloritos = np.zeros((nrgb,3))
    # coloritos[i_,0] = ...
    for i_ in range(0,nrgb):
      coloritos[i_,0] = r[i_]
      coloritos[i_,1] = g[i_]
      coloritos[i_,2] = b[i_]
    coloritomap = ListedColormap(coloritos)
    # ------------------------------------------------------------------------------
    return coloritomap
    # ------------------------------------------------------------------------------