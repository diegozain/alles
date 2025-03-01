import pyvista as pv
import numpy as np
import sys
# --------------------------------------------------------------------
def pointclic(mesh, picked_point):
  pointclicked = picked_point.GetPointId()
  print(f" ⦿ ⦿ .. {pointclicked:d}")
  
  global aii
  plotter.remove_actor(aii)
  
  # ▢
  aii=plotter.add_points(
  points[pointclicked,:],
  style="points_gaussian",
  # render_points_as_spheres=True,
  point_size=30,
  color="black"
  )

def closeclic(flag):
  print("⛔ bye bye")
  sys.exit(0)

instructxt = " click on a point and see its index in the console\n click the red square to exit"
# --------------------------------------------------------------------
# 10⁶ points makes the plot go 🐌 in my laptop
# 10⁵ points works 🚀 in my laptop
num_points = 10

np.random.seed(42)
points = np.random.random((num_points, 3))
# --------------------------------------------------------------------
point_cloud = pv.PolyData(points)
# --------------------------------------------------------------------
# add point indices as scalars
point_cloud["index"] = np.arange(num_points)
# initialize global variable to store clicked point index (??)
pointclicked = 0
# --------------------------------------------------------------------
# 📊
plotter = pv.Plotter()

# ▢
aii = plotter.add_points(
points[pointclicked,:],
style="points_gaussian",
# render_points_as_spheres=True,
point_size=30,
color="black"
)

# ⦿
plotter.enable_point_picking(
  callback=pointclic,
  show_message=False,
  use_picker=True,
)

# ····..
# ····..
plotter.add_mesh(
  point_cloud,
  render_points_as_spheres=True,
  point_size=30,
  scalars="index",
  cmap="viridis",
  pickable=True
)

# 📔
plotter.add_text(
  instructxt,
  font_size=12,
  position="upper_left"
)

# ⛔
plotter.add_checkbox_button_widget(
  callback=closeclic,
  value=False,
  color_on="white",
  color_off="red",
  background_color="black",
  size=40,
  position=(10.0,10.0)
)

plotter.show()
# --------------------------------------------------------------------
