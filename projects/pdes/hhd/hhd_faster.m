clear
close all
addpath('src/');
addpath('../../graph-alg/mesher/src/');
% ------------------------------------------------------------------------------
% diego domenzain
% 
% compute the helmholtz-hodge decomposition.
% 
% based on:
% The Natural Helmholtz-Hodge Decomposition for Open-Boundary Flow Analysis. 
% Harsh Bhatia, Valerio Pascucci, Peer-Timo Bremer. 
% IEEE Transactions on Visualization and Computer Graphics, 2014.
% ------------------------------------------------------------------------------
nx=50;
ny=50;
% ------------------------------------------------------------------------------
y=linspace(-1, 1,ny);
x=linspace(-1, 1,nx);

dx=x(2)-x(1);
dy=y(2)-y(1);
% ------------------------------------------------------------------------------
[xx,yy] = meshgrid(x,y);

% % stripes like the pics in readme
% ux =   sin(pi*xx).*cos(pi*yy) + sin(pi*yy).*cos(pi*xx) + 0.5*ones(ny,nx);
% uy = - sin(pi*xx).*cos(pi*yy) - sin(pi*yy).*cos(pi*xx) + 0.5*ones(ny,nx);
% 
% % stripes like the example in the paper (not referenced in the readme)
% ux =   sin(pi*xx).*cos(pi*yy) + sin(pi*yy).*cos(pi*xx) + 0.5*ones(ny,nx);
% uy = - sin(pi*yy).*cos(pi*xx) + sin(pi*xx).*cos(pi*yy) - ones(ny,nx);

% stripes like in the readme but no dc
ux =   sin(pi*xx).*cos(pi*yy) + sin(pi*yy).*cos(pi*xx);
uy = - sin(pi*xx).*cos(pi*yy) - sin(pi*yy).*cos(pi*xx);

clear xx yy;
% ------------------------------------------------------------------------------
figure;

subplot(1,2,1)
fancy_imagesc(ux,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('True ux')
simple_figure()

subplot(1,2,2)
fancy_imagesc(uy,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('True uy')
simple_figure()
% ------------------------------------------------------------------------------
% 
% 
%                       helmholtz-hodge decomposition
% 
% 
% ------------------------------------------------------------------------------
[xx,yy] = meshgrid(x,y);

L_=zeros(nx*ny,nx*ny);

for ix_=1:nx
 for iy_=1:ny
  % compute green's solution with source at (xo,yo)
  xo=x(ix_);
  yo=y(iy_);
  g = - (1/2/pi)*log( sqrt( ( xx-xo ).^2 + ( yy-yo ).^2 ) );
  % manage singularity at source location
  g(iy_,ix_) = 0;
  % get integral operator
  il = sub2ind([ny,nx],iy_,ix_);
  L_(:,il) = g(:);
 end
end
% ------------------------------------------------------------------------------
tic;
[phi_,psi_] = hhd_L_(ux,uy,dx,dy,nx,ny,L_);
% ------------------------------------------------------------------------------
% rebuild
phi_x = differentiate_plane(phi_.',dx);
phi_x = phi_x.';
phi_y = differentiate_plane(phi_,dy);

psi_x = differentiate_plane(psi_.',dx);
psi_x = psi_x.';
psi_y = differentiate_plane(psi_,dy);

hx = ux - phi_x - psi_y;
hy = uy - phi_y + psi_x;

ux = phi_x + psi_y + hx;
uy = phi_y - psi_x + hy;
% ------------------------------------------------------------------------------
% get derivatives 
hxx = differentiate_plane(hx.',dx);
hxx = hxx.';
hyy = differentiate_plane(hy,dy);

hxy = differentiate_plane(hx,dy);
hyx = differentiate_plane(hy.',dx);
hyx = hyx.';

divh = hxx + hyy;
roth = hyx - hxy;

toc;
% ------------------------------------------------------------------------------
% 
% 
%                                visualize
% 
% 
% ------------------------------------------------------------------------------
figure;

subplot(4,2,1)
fancy_imagesc(ux,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('ux')
simple_figure()

subplot(4,2,2)
fancy_imagesc(uy,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('uy')
simple_figure()

subplot(4,2,3)
fancy_imagesc(phi_x,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Curl-free x')
simple_figure()

subplot(4,2,4)
fancy_imagesc(phi_y,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Curl-free y')
simple_figure()

subplot(4,2,5)
fancy_imagesc(psi_y,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Div-free x')
simple_figure()

subplot(4,2,6)
fancy_imagesc(-psi_x,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Div-free y')
simple_figure()

subplot(4,2,7)
fancy_imagesc(hx,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Harmonic x')
simple_figure()

subplot(4,2,8)
fancy_imagesc(hy,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Harmonic y')
simple_figure()
% ------------------------------------------------------------------------------
figure;

subplot(2,4,1)
fancy_imagesc(ux(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,2)
fancy_imagesc(phi_x(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,3)
fancy_imagesc(psi_y(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,4)
fancy_imagesc(hx(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,5)
fancy_imagesc(uy(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,6)
fancy_imagesc(phi_y(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,7)
fancy_imagesc(-psi_x(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,8)
fancy_imagesc(hy(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()
% ------------------------------------------------------------------------------
figure;

subplot(1,2,1)
fancy_imagesc(phi_,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(1,2,2)
fancy_imagesc(psi_,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()
% ------------------------------------------------------------------------------
% % }
