function [phi_,psi_] = hhd_L_(ux,uy,dx,dy,nx,ny,L_)
% diego domenzain
% @ colorado school of mines, spring 2021
% 
% numerical implementation of the Helmholtz-Hodge decomposition.
% ------------------------------------------------------------------------------
% compute derivatives
uxx = differentiate_plane(ux.',dx);
uxx = uxx.';
uyy = differentiate_plane(uy,dy);

uxy = differentiate_plane(ux,dy);
uyx = differentiate_plane(uy.',dx);
uyx = uyx.';

divu= (uxx + uyy)*(dx*dy);
rotu= (uyx - uxy)*(dx*dy);

clear uxx uyy uxy uyx;
% ------------------------------------------------------------------------------
divu=divu(:);
rotu=rotu(:);
% ------------------------------------------------------------------------------
% solve
% ------------------------------------------------------------------------------
phi_ = L_*divu;
psi_ = L_*rotu;

phi_ = reshape(phi_,ny,nx);
psi_ = reshape(psi_,ny,nx);
end