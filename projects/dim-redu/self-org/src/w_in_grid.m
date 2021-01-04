function w = w_in_grid(v,u,ny,nx)
% initialize weights of a graph-grid (ny by nx)
% as a lattice in a plane spanned by 
% v and u.
% v and u can be in R^n.
% --
% w is a matrix of size (n by nx*ny),
% rows are entries of R^n,
% columns are nodes in the grid.
% ---------------
% center of grid
% nx_ = fix(nx*0.5);
% ny_ = fix(ny*0.5);
% ---------------
% coefficients for li. comb.
cx = linspace(-1,1,nx);
cy = linspace(-1,1,ny);
% ---------------
% create weights
% dimension of weights
nw = numel(u);
w = zeros(nw,nx*ny);
for i_=1:nx
  for j_=1:ny
    iw = sub2ind([ny,nx],j_,i_);
    w(:,iw) = cx(i_)*u + cy(j_)*v;
  end
end
end