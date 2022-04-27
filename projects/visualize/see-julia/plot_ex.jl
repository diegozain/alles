# ------------------------------------------------------------------------------
using Plots

x=1:10;
y=rand(10);
fig=plot(x,y);
display(fig);

function f(x,y)
r=sqrt(x^2 + y^2)
return cos(r)/(1+r)
end

x=range(0,2*pi,length=30);

fig=heatmap(x,x,f,c=:thermal);
display(fig);
# ------------------------------------------------------------------------------
#
#                      https://makie.juliaplots.org/stable/
#
#   https://makie.juliaplots.org/stable/documentation/colors/index.html#colors
#
# ------------------------------------------------------------------------------
# import Pkg; Pkg.add("GLMakie")
using GLMakie

r = LinRange(-1, 1, 100);
cube = [(x.^2 + y.^2 + z.^2) for x = r, y = r, z = r];
fig = contour(cube, alpha=0.5);
display(fig);

xs = [1, 2, 3, 1, 2, 3, 1, 2, 3];
ys = [1, 1, 1, 2, 2, 2, 3, 3, 3];
zs = [1, 2, 3, 4, 5, 6, 7, 8, NaN];
fig=heatmap(xs, ys, zs,colormap=:darktest);
display(fig);
# ------------------------------------------------------------------------------
# import Pkg; Pkg.add("CairoMakie")
using CairoMakie
using ElectronDisplay
ElectronDisplay.CONFIG.single_window = true;

x = range(0, 10, length=100);
y = sin.(x);
fig=lines(x, y);
display(fig);
# ------------------------------------------------------------------------------
