function [g_,g] = graph_torus(ny,nx)
% build incidence relations of an,
% ny by nx grid.
g = 1:(nx*ny);
g = reshape(g,ny,nx);
% ----------
% neighbors
% (counter clockwise like polar coord)
% ----------
g_right = [g(:,2:nx) , g(:,1)];
g_upper = [g(ny,:) ; g(1:(ny-1),:)];
g_left  = [g(:,nx) , g(:,1:(nx-1))];
g_down  = [g(2:ny,:) ; g(1,:)];
% this g is a cube (ny x nx x 5) where:
% the first (ny x nx) plane labels of the nodes,
% the second (ny x nx) plane labels the right neighbors,
% etc.
g = cat(3,g,g_right,g_upper,g_left,g_down);
clear g_right g_upper g_left g_down;
% -----------
% make neighbors list
% -----------
% this g_ has the same info as g above but now as 
% a (nx*ny x 4) matrix where each row corresponds to a node,
% and each column to a neighbor, 
% i.e. g_(i,:) are all the neighbors of node i.
g_ = zeros(nx*ny,4);
for i_=1:(nx*ny) % this for loop is parforable.
  [iy,ix] = ind2sub([ny,nx],i_);
  g_(i_,:) = [g(iy,ix,2) , g(iy,ix,3) , g(iy,ix,4) , g(iy,ix,5)];
end
% all of g and g_ are indicies:
g = uint32(g);
g_ = uint32(g_);
end