function p = ray_path(r,s,grad_u,h,nx)
%
% ray r ------> s
%
% r and s are in index notation,
% i.e. an integer corresponding to an entry in the xz-plane.

% bring to matrix notation
r = p2ij(r,nx);
s = p2ij(s,nx);

% for gradient descent
tau = h;

p = r;
% loop
while norm(p(end,:)-s) > 2e+0
% for i=1:1000
  % get (x,z) coordinate of point in binned indicies
  % p_binned = [ binning(x,p(end,1)) , binning(z,p(end,2)) ];
  p_binned = fix(p(end,:));

  % evaluate gradient at point
  p_ = ij2p( p_binned , nx);
  g = grad_u(p_,:);

  % update
  p_ = p(end,:) - tau*g;
  p = [p ; p_]; 
end

end