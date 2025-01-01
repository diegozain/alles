# 🎨🐍

plot in python.

```python
# 🟪🟧🟥🟨⚫
from coloritobar import coloritobar
colorito = coloritobar()

vmin = matrix.min()
vmid = 0.0
vmax = matrix.max()
coloritomap = colorito.build(vmin,vmid,vmax,800)

fig, ax = plt.subplots()
im=ax.imshow(matrix, aspect='auto',cmap=coloritomap,vmin=vmin,vmax=vmax)
fig.colorbar(im,ax=ax,orientation="horizontal",label="Colormap",pad=0.2,fraction=0.05)
plt.show()
```