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

# ------------------------------------------------------
# import Pkg; Pkg.add("GLMakie")

using GLMakie
r = LinRange(-1, 1, 100);
cube = [(x.^2 + y.^2 + z.^2) for x = r, y = r, z = r];
fig = contour(cube, alpha=0.5);
display(fig);

# ------------------------------------------------------
# import Pkg; Pkg.add("CairoMakie")

using CairoMakie

x = range(0, 10, length=100);
y = sin.(x);
fig=lines(x, y);

# ------------------------------------------------------
