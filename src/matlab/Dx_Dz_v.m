function [Dx,Dz] = Dx_Dz_v(n,m)
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
  fwd = zeros(n,n*m);
  dz_fwd = sparse(fwd);
  
  % center
  ctd = [-ones(1,n*m-n) ones(1,n*m-n)];
  I = 1:n*m-n;
  I = repmat(I,1,2);
  J = [1:(n*m-n) ((n+1):n*m)];
  dz_ctd = sparse(I,J,ctd);
  
  % % center
  % ctd = [-ones(1,n*m-2*n) ones(1,n*m-2*n) zeros(1,n*m-2*n)];
  % I = 1:n*m-2*n;
  % I = repmat(I,1,3);
  % J = [1:(n*m-2*n) (n+1):(n*m-2*n+n) (n+n+1):(n*m-2*n+n+n)];
  % dz_ctd = sparse(I,J,ctd);
  % % bwd
  % bwd = zeros(n,n*m);
  % dz_bwd = sparse(bwd);
  
  % together
  Dz = [dz_fwd; dz_ctd];
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
  
  fwd = zeros(1,n);  
  dx_fwd = sparse(fwd);
  
  ctd = [-ones(1,n-1) ones(1,n-1)];
  I = 1:n-1;
  I = repmat(I,1,2);
  J = [(1:n-1) (1+1:n-1+1)];
  dx_ctd = sparse(I,J,ctd);
  
  % ctd = [-ones(1,n-2) ones(1,n-2) zeros(1,n-2)];
  % I = 1:n-2;
  % I = repmat(I,1,3);
  % J = [(1:n-2) (1+1:n-2+1) (1+1+1:n-2+1+1)];
  % dx_ctd = sparse(I,J,ctd);
  % %
  % bwd = zeros(1,n);
  % dx_bwd = sparse(bwd);
  
  Dx_ = [dx_fwd; dx_ctd];
end