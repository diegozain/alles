function [phi_,psi_] = hhd(ux,uy,x,y)

[ny,nx] = size(ux);
[xx,yy] = meshgrid(x,y);
dx = x(2)-x(1);
dy = y(2)-y(1);
% ------------------------------------------------------------------------------

% compute derivatives
uxx = differentiate_plane(ux.',dx);
uxx = uxx.';
uyy = differentiate_plane(uy,dy);

uxy = differentiate_plane(ux,dy);
uyx = differentiate_plane(uy.',dx);
uyx = uyx.';

divu = uxx + uyy;
rotu = uyx - uxy;

clear uxx, uyy, uxy, uyx
% ------------------------------------------------------------------------------
% solve
% ------------------------------------------------------------------------------

phi_=zeros(ny,nx);
psi_=zeros(ny,nx);

for ix_=1:nx
 for iy_=1:ny
  % compute green's solution with source at (xo,yo)
  xo=x(ix_);
  yo=y(iy_);
  g = - (1/2/pi)*log( sqrt( ( xx-xo ).^2 + ( yy-yo ).^2 ) );
  % manage singularity at source location
  g(iy_,ix_) = 0;
  % get phi_ values (integral is clunky: just sum)
  a = g.*divu;
  phi_(iy_,ix_) = sum(a(:))*dx*dy;
  % get psi_ values (integral is clunky: just sum)
  b = g.*rotu;
  psi_(iy_,ix_) = sum(b(:))*dx*dy;
 end
end

end