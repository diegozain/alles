import pyvista as pv
import numpy as np
import sys
# --------------------------------------------------------------------
# Define the callback function to handle picking events
def pointclic(mesh, picked_point):
  print(f" ⦿ ⦿ .. {picked_point.GetPointId()}")

def closeclic(flag):
  print("⛔ bye bye")
  sys.exit(0)

textplot = " click on a point and\n see its index in the console.\n click the red button to exit."
# --------------------------------------------------------------------
# 10⁶ points makes the plot go 🐌 in my laptop
# 10⁵ points works 🚀 in my laptop
num_points = 500000

np.random.seed(42)
points = np.random.random((num_points, 3))
# --------------------------------------------------------------------
point_cloud = pv.PolyData(points)
# --------------------------------------------------------------------
# add point indices as scalars for easy identification (??)
point_cloud["index"] = np.arange(num_points)
# Initialize global variable to store clicked point index
selected_index = None
# --------------------------------------------------------------------
# 📊
plotter = pv.Plotter()

plotter.add_text(
  textplot,
  font_size=12,
  position="upper_left"
)

plotter.add_mesh(
  point_cloud,
  render_points_as_spheres=True,
  point_size=30,
  scalars="index",
  cmap="viridis",
  pickable=True
)

# ⦿
plotter.enable_point_picking(
  callback=pointclic,
  show_message=False,
  use_picker=True,
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

# plotter
# ['_BasePlotter__before_close_callback', '__class__', '__del__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_add_legend_label', '_before_close_callback', '_check_has_ren_win', '_check_rendered', '_clear_picking_representations', '_clear_ren_win', '_closed', '_datasets', '_first_time', '_has_background_layer', '_id_name', '_image_depth_null', '_image_scale', '_init_click_picking_callback', '_initialized', '_make_render_window_current', '_on_first_render_request', '_on_render_callbacks', '_picked_actor', '_picked_block_index', '_picked_cell', '_picked_mesh', '_picked_point', '_picker_in_use', '_picking_left_clicking_observer', '_picking_right_clicking_observer', '_picking_text', '_prep_for_close', '_rendered', '_save_image', '_scalar_bar_slot_lookup', '_scalar_bar_slots', '_scalar_bars', '_suppress_rendering', '_theme', '_validate_picker_not_in_use', '_window_size_unset', 'actors', 'add_actor', 'add_affine_transform_widget', 'add_arrows', 'add_axes', 'add_axes_at_origin', 'add_background_image', 'add_blurring', 'add_bounding_box', 'add_box_axes', 'add_box_widget', 'add_camera3d_widget', 'add_camera_orientation_widget', 'add_chart', 'add_checkbox_button_widget', 'add_composite', 'add_cursor', 'add_floor', 'add_key_event', 'add_legend', 'add_legend_scale', 'add_light', 'add_line_widget', 'add_lines', 'add_logo_widget', 'add_measurement_widget', 'add_mesh', 'add_mesh_clip_box', 'add_mesh_clip_plane', 'add_mesh_isovalue', 'add_mesh_slice', 'add_mesh_slice_orthogonal', 'add_mesh_slice_spline', 'add_mesh_threshold', 'add_north_arrow_widget', 'add_on_render_callback', 'add_orientation_widget', 'add_plane_widget', 'add_point_labels', 'add_point_scalar_labels', 'add_points', 'add_ruler', 'add_scalar_bar', 'add_silhouette', 'add_slider_widget', 'add_sphere_widget', 'add_spline_widget', 'add_text', 'add_text_slider_widget', 'add_timer_event', 'add_title', 'add_volume', 'add_volume_clip_plane', 'background_color', 'bounds', 'box_clipped_meshes', 'box_widgets', 'button_widgets', 'camera', 'camera3d_widgets', 'camera_position', 'camera_set', 'camera_widgets', 'center', 'clear', 'clear_actors', 'clear_box_widgets', 'clear_button_widgets', 'clear_camera3d_widgets', 'clear_camera_widgets', 'clear_events_for_key', 'clear_line_widgets', 'clear_logo_widgets', 'clear_measure_widgets', 'clear_on_render_callbacks', 'clear_plane_widgets', 'clear_slider_widgets', 'clear_sphere_widgets', 'clear_spline_widgets', 'click_position', 'close', 'deep_clean', 'disable', 'disable_3_lights', 'disable_anti_aliasing', 'disable_depth_of_field', 'disable_depth_peeling', 'disable_eye_dome_lighting', 'disable_hidden_line_removal', 'disable_parallel_projection', 'disable_picking', 'disable_shadows', 'disable_ssao', 'disable_stereo_render', 'distance_widgets', 'enable', 'enable_2d_style', 'enable_3_lights', 'enable_anti_aliasing', 'enable_block_picking', 'enable_cell_picking', 'enable_custom_trackball_style', 'enable_depth_of_field', 'enable_depth_peeling', 'enable_element_picking', 'enable_eye_dome_lighting', 'enable_fly_to_right_click', 'enable_geodesic_picking', 'enable_hidden_line_removal', 'enable_horizon_picking', 'enable_image_style', 'enable_joystick_actor_style', 'enable_joystick_style', 'enable_lightkit', 'enable_mesh_picking', 'enable_parallel_projection', 'enable_path_picking', 'enable_point_picking', 'enable_rectangle_picking', 'enable_rectangle_through_picking', 'enable_rectangle_visible_picking', 'enable_rubber_band_2d_style', 'enable_rubber_band_style', 'enable_shadows', 'enable_ssao', 'enable_stereo_render', 'enable_surface_point_picking', 'enable_terrain_style', 'enable_trackball_actor_style', 'enable_trackball_style', 'enable_zoom_style', 'export_gltf', 'export_html', 'export_obj', 'export_vrml', 'export_vtksz', 'fly_to', 'fly_to_mouse_position', 'generate_orbital_path', 'get_default_cam_pos', 'get_image_depth', 'get_pick_position', 'hide_axes', 'hide_axes_all', 'image', 'image_depth', 'image_scale', 'image_scale_context', 'image_transparent_background', 'import_3ds', 'import_gltf', 'import_obj', 'import_vrml', 'increment_point_size_and_line_width', 'iren', 'isometric_view', 'isometric_view_interactive', 'isovalue_meshes', 'key_press_event', 'last_image', 'last_image_depth', 'last_update_time', 'last_vtksz', 'left_button_down', 'legend', 'length', 'line_widgets', 'link_views', 'logo_widgets', 'mesh', 'meshes', 'mouse_position', 'notebook', 'off_screen', 'open_gif', 'open_movie', 'orbit_on_path', 'parallel_projection', 'parallel_scale', 'pick_click_position', 'pick_mouse_position', 'pickable_actors', 'picked_actor', 'picked_block_index', 'picked_cell', 'picked_cells', 'picked_geodesic', 'picked_horizon', 'picked_mesh', 'picked_path', 'picked_point', 'plane_clipped_meshes', 'plane_sliced_meshes', 'plane_widgets', 'remove_actor', 'remove_all_lights', 'remove_background_image', 'remove_blurring', 'remove_bounding_box', 'remove_bounds_axes', 'remove_chart', 'remove_environment_texture', 'remove_floors', 'remove_legend', 'remove_scalar_bar', 'ren_win', 'render', 'render_window', 'renderer', 'renderers', 'reset_camera', 'reset_camera_clipping_range', 'reset_key_events', 'save_graphic', 'scalar_bar', 'scalar_bars', 'scale', 'screenshot', 'set_background', 'set_chart_interaction', 'set_color_cycler', 'set_environment_texture', 'set_focus', 'set_position', 'set_scale', 'set_viewup', 'shape', 'show', 'show_axes', 'show_axes_all', 'show_bounds', 'show_grid', 'slider_widgets', 'sphere_widgets', 'spline_sliced_meshes', 'spline_widgets', 'store_click_position', 'store_mouse_position', 'subplot', 'suppress_rendering', 'theme', 'threshold_meshes', 'title', 'track_click_position', 'track_mouse_position', 'unlink_views', 'untrack_click_position', 'untrack_mouse_position', 'update', 'update_bounds_axes', 'update_coordinates', 'update_scalar_bar_range', 'update_scalars', 'view_isometric', 'view_vector', 'view_xy', 'view_xz', 'view_yx', 'view_yz', 'view_zx', 'view_zy', 'where_is', 'window_size', 'window_size_context', 'write_frame', 'zoom_camera']

# picked_point
# ['AddObserver', 'AddPickList', 'BreakOnError', 'DebugOff', 'DebugOn', 'DeletePickList', 'FastDelete', 'GetActor', 'GetActor2D', 'GetActors', 'GetAddressAsString', 'GetAssembly', 'GetClassName', 'GetCommand', 'GetCompositeDataSet', 'GetDataSet', 'GetDebug', 'GetFlatBlockIndex', 'GetGlobalWarningDisplay', 'GetIsInMemkind', 'GetMTime', 'GetMapper', 'GetMapperPosition', 'GetNumberOfGenerationsFromBase', 'GetNumberOfGenerationsFromBaseType', 'GetObjectDescription', 'GetObjectName', 'GetPath', 'GetPickFromList', 'GetPickList', 'GetPickPosition', 'GetPickedPositions', 'GetPointId', 'GetProp3D', 'GetProp3Ds', 'GetPropAssembly', 'GetReferenceCount', 'GetRenderer', 'GetSelectionPoint', 'GetTolerance', 'GetUseCells', 'GetUsingMemkind', 'GetViewProp', 'GetVolume', 'GlobalWarningDisplayOff', 'GlobalWarningDisplayOn', 'HasObserver', 'InitializeObjectBase', 'InitializePickList', 'InvokeEvent', 'IsA', 'IsTypeOf', 'Modified', 'NewInstance', 'Pick', 'Pick3DPoint', 'Pick3DRay', 'PickFromListOff', 'PickFromListOn', 'Register', 'RemoveAllObservers', 'RemoveObserver', 'RemoveObservers', 'SafeDownCast', 'SetDebug', 'SetGlobalWarningDisplay', 'SetMemkindDirectory', 'SetObjectName', 'SetPath', 'SetPickFromList', 'SetReferenceCount', 'SetTolerance', 'SetUseCells', 'UnRegister', 'UseCellsOff', 'UseCellsOn', 'UsesGarbageCollector', '__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__this__', '__vtkname__', 'override']