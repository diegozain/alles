function [Dx,Dz] = Dx_Dz_s(n,m)
  % diego domenzain
  % fall 2021 @ CSM
  % ----------------------------------------------------------------------------
  % 
  % ----------------------------------------------------------------------------
  Dx = DX(n,m);
  Dz = DZ(n,m);
end

% --------------------------------------------
%   z derivative
% --------------------------------------------

function Dz = DZ(n,m)
  
  % fwd
  fwd = [-ones(1,n*m-1*n) ones(1,n*m-1*n)];
  I = 1:n*m-1*n;
  I = repmat(I,1,2);
  J = [(1:n*m-1*n) (n+1:n*m-1*n+n)];
  dz_fwd = sparse(I,J,fwd);
  
  % bwd
  bwd = zeros(n,n*m);
  dz_bwd = sparse(bwd);
  
  % together
  Dz = [dz_fwd; dz_bwd];
end

% --------------------------------------------
%   x derivative
% --------------------------------------------

% first degree derivative, 
% second degree accurate.
% 
% fwd:             −1.5 	 2 	−0.5
% ctd:        −0.5 	 0 	 0.5
% bwd:   0.5 	−2 	 1.5

function Dx = DX(n,m)
  % takes first derivative on n points,
  % and does this for m times,
  %
  %                 n
  %        .-------------------.
  %        |                   |
  %        |                   |
  %    m   |       DX          |
  %        |                   |
  %        |                   |
  %        .-------------------.
  %
  %
  %       .-----.
  %       |     |
  %       |     | n     f is (n*m x 1)
  %       |     | .
  %       |  f  | .
  %       |     | .
  %       |     | n
  %       |     |
  %       .-----.
  %
  % DX * f = dx(f)
  %
  % where dx is in the direction of n.
  
  Dx = D_x(n);
  id = speye(m);
  
  Dx = kron(id,Dx);
end

function Dx_ = D_x(n)
  %
  % second degree accurate,
  % first derivative matrix for
  % 1D discretization on
  % n pts.
  
  fwd = [-ones(1,n-1) ones(1,n-1)];
  I = 1:n-1;
  I = repmat(I,1,2);
  J = [(1:n-1) (1+1:n-1+1)];
  dx_fwd = sparse(I,J,fwd);
  
  bwd = zeros(1,n);
  dx_bwd = sparse(bwd);
  
  Dx_ = [dx_fwd; dx_bwd];
end