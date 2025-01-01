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

  def build(self, minval, midval, maxval, nrgb):
    # ------------------------------------------------------------------------------
    #
    #                                    🔴🔶🔷
    #                                 build colormap
    #                                       σ
    # 
    # h="B4FBB8"
    # a=tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
    # ------------------------------------------------------------------------------
    # 🟣🔴 214
    coloritosr = np.asarray([9/255,160/255,0.9,113/255,176/255])
    coloritosg = np.asarray([9/255,136/255,0.9,133/255,90/255])
    coloritosb = np.asarray([11/255,17/255,0.9,145/255,58/255])
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