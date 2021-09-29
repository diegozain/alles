function dij = deltas3d(ig,neigh_mesh,iyxz,iyxz_,x,y,z)
% diego domenzain
% sep 2021
% ------------------------------------------------------------------------------
%
%      δij = Δij / Δ⟂ij
%
%          .|-----.
%         / |    /|←Δij1
%        /  |   / |             Δij = Δij1 · Δij2
%       .------.  |↓Δ⟂ij
%   ----|---•--| ----• iyxz_
%       |   |  | /
%       .------./←Δij2
%           |
% ------------------------------------------------------------------------------
% neighbor types. these are in the 4rth column of robin_xyz.
%      2  6
%      | /
% 3 -- i -- 1
%    / |
%   5  4
%
%  1 2 3 4 5  6
%  → ↑ ← ↓ ⦿ ⊛
%
% the numbers 1, 2, 3, 4, 5 and 6,
% represent neighbors,
% right, up, left, down, front and back.
% ------------------------------------------------------------------------------
% get the type of neighobr iyxz_ is to iyxz
inei =0;
click=1;
while (inei<7 && click==1)
  inei  = inei+1;
  iyxz__= neigh_mesh(ig,inei);
  if (iyxz__==iyxz_)
    click = 0;
  end
end
% ------------------------------------------------------------------------------
% 🍹 trick!
% there is a problem accessing Δij1, Δij2, and Δ⟂ij at the boundaries,
% so here we fix that by adding elements to x, y, and z:
%
%     1        ix            nx
%     .---------.-------------.
%     2       ix+1          nx+1
%  .--.---------.-------------.--.      ← added elements
%  1                           nx+2
nx=numel(x);
ny=numel(y);
nz=numel(z);

[ix,iy,iz]    = get_ixyz(iyxz,nx,ny,nz);
[ix_,iy_,iz_] = get_ixyz(iyxz_,nx,ny,nz);

ix =ix+1;
iy =iy+1;
iz =iz+1;

ix_=ix_+1;
iy_=iy_+1;
iz_=iz_+1;

x = [x(1) - (x(2)-x(1)) ; x ; x(nx) + (x(nx)-x(nx-1))];
y = [y(1) - (y(2)-y(1)) ; y ; y(ny) + (y(ny)-y(ny-1))];
z = [z(1) - (z(2)-z(1)) ; z ; z(nz) + (z(nz)-z(nz-1))];
% ------------------------------------------------------------------------------
%  1 2 3 4 5  6
%  → ↑ ← ↓ ⦿ ⊛

%  →
if (inei==1)
  dij_ = x(ix+1) - x(ix);
  dij1 = (y(iy+1)- y(iy-1))*0.5;
  dij2 = (z(iz+1)- z(iz-1))*0.5;
%  ↑
elseif (inei==2)
  dij_ = z(iz-1)  - z(iz);
  dij1 = (y(iy+1) - y(iy-1))*0.5;
  dij2 = (x(ix+1) - x(ix-1))*0.5;
%  ←
elseif (inei==3)
  dij_ = x(ix-1) - x(ix);
  dij1 = (y(iy+1) - y(iy-1))*0.5;
  dij2 = (z(iz+1) - z(iz-1))*0.5;
%  ↓
elseif (inei==4)
  dij_ = z(iz+1)  - z(iz);
  dij1 = (y(iy+1) - y(iy-1))*0.5;
  dij2 = (x(ix+1) - x(ix-1))*0.5;
%  ⦿
elseif (inei==5)
  dij_ = y(iy+1)  - y(iy);
  dij1 = (z(iz+1) - z(iz-1))*0.5;
  dij2 = (x(ix+1) - x(ix-1))*0.5;
%  ⊛
elseif (inei==6)
  dij_ = y(iy)  - y(iy-1);
  dij1 = (z(iz+1) - z(iz-1))*0.5;
  dij2 = (x(ix+1) - x(ix-1))*0.5;
end

dij = abs((dij1*dij2) / dij_);
end
