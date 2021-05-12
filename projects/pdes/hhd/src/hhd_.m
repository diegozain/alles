function [phi_,psi_] = hhd_(ux,uy,x,y)
% diego domenzain
% @ colorado school of mines, spring 2021
% 
% numerical implementation of the Helmholtz-Hodge decomposition.
% loosely based on:
% The Natural Helmholtz-Hodge Decomposition for Open-Boundary Flow Analysis. 
% Harsh Bhatia, Valerio Pascucci, Peer-Timo Bremer. 
% IEEE Transactions on Visualization and Computer Graphics, 2014.
% ------------------------------------------------------------------------------
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

clear uxx uyy uxy uyx;
% ------------------------------------------------------------------------------
% solve
% ------------------------------------------------------------------------------
% uses Fortran-Mex code:
[phi_,psi_] = hhd_solve(x,y,divu,rotu);
end