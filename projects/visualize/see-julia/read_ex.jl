using GLMakie
using ColorSchemes
# https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/

# using CairoMakie
# using ElectronDisplay
# ElectronDisplay.CONFIG.single_window = true;
# ------------------------------------------------------------------------------
if Sys.iswindows()
 include("..\\..\\..\\src\\julia\\readwrite\\read_binf64.jl")
 include("..\\..\\..\\src\\julia\\readwrite\\read_bini32.jl")

 include("..\\..\\..\\src\\julia\\mesher3d\\get_ixyz.jl")
 include("..\\..\\..\\src\\julia\\mesher3d\\n_g2m_3d_.jl")
 include("..\\..\\..\\src\\julia\\mesher3d\\g2m_m2g_3d.jl")
 include("..\\..\\..\\src\\julia\\mesher3d\\cubify_graph.jl")
elseif Sys.islinux()
 include("../../../src/julia/readwrite/read_binf64.jl")
 include("../../../src/julia/readwrite/read_bini32.jl")

 include("../../../src/julia/mesher3d/get_ixyz.jl")
 include("../../../src/julia/mesher3d/n_g2m_3d_.jl")
 include("../../../src/julia/mesher3d/g2m_m2g_3d.jl")
 include("../../../src/julia/mesher3d/cubify_graph.jl")
end
# ------------------------------------------------------------------------------
#
#
#                                  📥📥📥
#
#
# ------------------------------------------------------------------------------
if Sys.iswindows()
 path_read="..\\..\\..\\..\\dcip\\field\\ejls\\jan22\\e6\\bin\\";
 path_read_="..\\..\\..\\..\\dcip\\field\\ejls\\mar22\\e20\\bin\\";
elseif Sys.islinux()
 path_read="../../../../dcip/field/ejls/jan22/e6/bin/";
 path_read_="../../../../dcip/field/ejls/mar22/e20/bin/";
end
# ------------------------------------------------------------------------------
nx= read_bini32(string(path_read,"x_size"),[3,1]);
nx=nx[1];
x = read_binf64(string(path_read,"x"),[nx,1]);
ny= read_bini32(string(path_read,"y_size"),[3,1]);
ny=ny[1];
y = read_binf64(string(path_read,"y"),[ny,1]);
nz= read_bini32(string(path_read,"z_size"),[3,1]);
nz=nz[1];
z = read_binf64(string(path_read,"z"),[nz,1]);

nmask3d= read_bini32(string(path_read,"mask3d_size"),[3,1]);
nmask3d=nmask3d[1];
mask3d = read_bini32(string(path_read,"mask3d"),[nmask3d,1]);

nmn= read_bini32(string(path_read,"mn_size"),[3,1]);
mn = read_bini32(string(path_read,"mn"),[nmn[1],nmn[2]]);
# ------------------------------------------------------------------------------
#                                   👷👷👷
# ------------------------------------------------------------------------------
n_g2m = n_g2m_3d_(mask3d,nx,ny,nz);
graph2mesh, mesh2graph = g2m_m2g_3d(mask3d,nx,ny,nz,n_g2m);
# ------------------------------------------------------------------------------
sigm_reco = read_binf64(string(path_read,"sigm_reco"),[n_g2m,1]);
sigm_reco_= read_binf64(string(path_read_,"sigm_reco"),[n_g2m,1]);
rho_reco = 1 ./ sigm_reco;
quot_reco = sigm_reco_ ./ sigm_reco;
# ------------------------------------------------------------------------------
sigm_reco3d_ = cubify_graph(sigm_reco,graph2mesh,nx,ny,nz,n_g2m);
rho_reco3d_ = cubify_graph(rho_reco,graph2mesh,nx,ny,nz,n_g2m);
quot_reco3d_ = cubify_graph(quot_reco,graph2mesh,nx,ny,nz,n_g2m);

sigm_reco2d_= sigm_reco3d_[25,:,:];
rho_reco2d_ = rho_reco3d_[25,:,:];
quot_reco2d_= quot_reco3d_[25,:,:];


x = dropdims(x;dims=2);
y = dropdims(y;dims=2);
z = dropdims(z;dims=2);
# ------------------------------------------------------------------------------
#
#
#                                   🎨🎨🎨
#
#
# ------------------------------------------------------------------------------
# fig=heatmap(x,z,sigm_reco2d_,colormap=:Egypt);
# display(fig);
# # ------------------------------------------------------------------------------
# fig=heatmap(x,z,rho_reco2d_,colormap=:Egypt);
# display(fig);
# # ------------------------------------------------------------------------------
# fig = Figure(resolution = (1200, 800))
# ax = Axis(fig[1, 1]; xlabel = "Length (m)", ylabel = "Depth (m)",title="÷");
# ax.yreversed = true;
# hmap=heatmap!(x,z,quot_reco2d_,colormap=:berlin);
# colsize!(fig.layout, 1, Aspect(1, 1.0));
# Colorbar(fig[1, 2], hmap; label = "( -- )", width = 15, ticksize = -10);
# colgap!(fig.layout, 10);
# display(fig);
# # ------------------------------------------------------------------------------
sigm_reco3d_=sigm_reco3d_[:,:,end:-1:1];
fig = Figure(resolution = (1200, 800))
ax = Axis3(fig[1, 1]; xlabel = "Width (m)", ylabel = "Length (m)",zlabel="Depth (m)",title="σ",aspect=:data);
# vol = volume!(y,x,z,sigm_reco3d_, colormap = :Egypt, transparency = true);
vol = volume!(y,x,z,sigm_reco3d_, colormap = :Egypt, algorithm = :iso, isorange = 0.05, isovalue = 0.01);
Colorbar(fig[1, 2], vol; label = "( S/m )", width = 15, ticksize = -10);
colgap!(fig.layout, -100);
display(fig);
# ------------------------------------------------------------------------------
# fig = volume(y,x,z,rho_reco3d_, colormap = :Egypt, transparency = true, figure = (; resolution = (1200, 800)), axis = (; type = Axis3, perspectiveness = 0.5, azimuth = 2.19, elevation = 0.57, aspect = (1, 1, 1)))
# ------------------------------------------------------------------------------
# fig = volume(y,x,z,quot_reco3d_, colormap = :Paired_10, transparency = true, figure = (; resolution = (1200, 800)), axis = (; type = Axis3, perspectiveness = 0.5, azimuth = 2.19, elevation = 0.57, aspect = (1, 1, 1)))
# # ------------------------------------------------------------------------------
# quot_reco3d_=quot_reco3d_[:,:,end:-1:1];
# fig = Figure(resolution = (1200, 800))
# ax = Axis3(fig[1, 1]; xlabel = "Width (m)", ylabel = "Length (m)",zlabel="Depth (m)",title="÷",aspect=:data);
# vol = volume!(y,x,z,quot_reco3d_, colormap = :Paired_10, algorithm = :iso, isorange = 0.001, isovalue = 0.5);
# Colorbar(fig[1, 2], vol; label = "( -- )", width = 15, ticksize = -10);
# colgap!(fig.layout, -100);
# ------------------------------------------------------------------------------
