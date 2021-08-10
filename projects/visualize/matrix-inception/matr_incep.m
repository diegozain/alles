clc
clear
close all
% ------------------------------------------------------------------------------
na = 10;
nb = 7;
nc = 13;
% ------------------------------------------------------------------------------
%        _____________
%      |  1           |
%      |  2           |
%  nb  |  3           |
%      |  etc         |
%      |______________|
%             na
% ------------------------------------------------------------------------------
% one for loop for two levels
mat2d = zeros(nb,na);
for iab = 1:na*nb
  % get b coordinate
  ib = mod(iab,nb);
  if (ib==0)
    ib=nb;
  end
  % iab = (ia - 1)*nb + ib
  % get a coordinate
  ia = ((iab-ib)/nb) + 1;

  mat2d(ib,ia) = iab;
end
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(mat2d)
colormap(rainbow2_cb(1))
xlabel('a')
ylabel('b')
title('two dimensional array')
simple_figure()
% ------------------------------------------------------------------------------
%             _____________
%           /             /|
%         /             /  |
%       /_____________/    |
%      |  1           |    |
%      |  2           |    |
%  nc  |  3           |    /
%      |  etc         |  / na
%      |______________|/
%            nb
% ------------------------------------------------------------------------------
% one for loop for two levels
mat3d = zeros(na,nb,nc);
mat3d_= zeros(na*nb*nc,4);
for iabc = 1:na*nb*nc
  % get c coordinate
  ic = mod(iabc,nc);
  if (ic==0)
    ic=nc;
  end
  % iabc = (ia-1)*nb*nc + (ib-1)*nc + ic  ... (*)
  ibc = mod(iabc,nb*nc);
  if (ibc==0)
    ibc=nb*nc;
  end
  % ibc = (ib-1)*nc + ic from (*)
  % get b coordinate
  ib = ((ibc-ic)/nc) + 1;
  % get a coordinate from (*)
  ia = ((iabc-ibc)/(nb*nc)) + 1;

  mat3d(ia,ib,ic) = iabc;
  mat3d_(iabc,:)  = [ia,ib,ic,iabc];
end
% ------------------------------------------------------------------------------
figure;
subplot(1,2,1)
fancy_imagesc(squeeze(mat3d(1,:,:)).')
colormap(rainbow2_cb(1))
caxis([1,na*nb*nc])
xlabel('b')
ylabel('c')
title('slice of 3d array')
simple_figure()

subplot(1,2,2)
scatter3(mat3d_(:,2),mat3d_(:,1),mat3d_(:,3),80*ones(na*nb*nc,1),mat3d_(:,4),'filled')
set(gca,'ZDir','reverse');
axis tight
colormap(rainbow2_cb(1));
caxis([1,na*nb*nc])
colorbar;
xlabel('b')
ylabel('a')
zlabel('c')
title('cube of 3d array')
simple_figure()
% ------------------------------------------------------------------------------
%             _____________
%           /             /|
%         /             /  |
%       /_____________/    |
%      |  1           |    |
%      |  2           |    |
%  nz  |  3           |    /
%      |  etc         |  / ny
%      |______________|/
%            nx
% ------------------------------------------------------------------------------
% nc --> nz
% nb --> nx
% na --> ny
% ------------------------------------------------------------------------------
nx=nb;
ny=na;
nz=nc;
% ------------------------------------------------------------------------------
% one for loop for two levels
mat3d = zeros(ny,nx,nz);
mat3d_= zeros(ny*nx*nz,4);
for iyxz = 1:ny*nx*nz
  % get z coordinate
  iz = mod(iyxz,nz);
  if (iz==0)
    iz=nz;
  end
  % iyxz = (iy-1)*nx*nz + (ix-1)*nz + iz  ... (*)
  ixz = mod(iyxz,nx*nz);
  if (ixz==0)
    ixz=nx*nz;
  end
  % ixz = (ix-1)*nz + iz from (*)
  % get x coordinate
  ix = ((ixz-iz)/nz) + 1;
  % get y coordinate from (*)
  iy = ((iyxz-ixz)/(nx*nz)) + 1;

  mat3d(iy,ix,iz) = iyxz;
  mat3d_(iyxz,:)  = [iy,ix,iz,iyxz];
end
% ------------------------------------------------------------------------------
figure;
subplot(1,2,1)
fancy_imagesc(squeeze(mat3d(1,:,:)).')
colormap(rainbow2_cb(1))
caxis([1,nx*ny*nz])
xlabel('x')
ylabel('z')
title('slice of 3d array')
simple_figure()

subplot(1,2,2)
scatter3(mat3d_(:,2),mat3d_(:,1),mat3d_(:,3),200*ones(ny*nx*nz,1),mat3d_(:,4),'filled')
set(gca,'ZDir','reverse');
axis tight;
axis image;
colormap(rainbow2_cb(1));
caxis([1,nx*ny*nz])
colorbar;
xlabel('x')
ylabel('y')
zlabel('z')
title('cube of 3d array')
simple_figure()
% ------------------------------------------------------------------------------
%             _____________
%           /             /|
%         /             /  |
%       /_____________/    |
%      |  1           |    |
%      |  2           |    |
%  nz  |  3    🎨    |    /
%      |  etc         |  / ny
%      |______________|/
%            nx
% ------------------------------------------------------------------------------
% nc --> nz
% nb --> nx
% na --> ny
% ------------------------------------------------------------------------------
x=linspace(0,1,13);
y=linspace(0,1,13);
z=linspace(0,1,13);

nx = numel(x);
ny = numel(y);
nz = numel(z);
% ------------------------------------------------------------------------------
% one for loop for two levels
mat3d = zeros(ny,nx,nz);
mat3d_= zeros(ny*nx*nz,4);
for iyxz = 1:ny*nx*nz
  % get z coordinate
  iz = mod(iyxz,nz);
  if (iz==0)
    iz=nz;
  end
  % iyxz = (iy-1)*nx*nz + (ix-1)*nz + iz  ... (*)
  ixz = mod(iyxz,nx*nz);
  if (ixz==0)
    ixz=nx*nz;
  end
  % ixz = (ix-1)*nz + iz from (*)
  % get x coordinate
  ix = ((ixz-iz)/nz) + 1;
  % get y coordinate from (*)
  iy = ((iyxz-ixz)/(nx*nz)) + 1;

  % val = exp(-(sum(sqrt( (x(ix)-0.5)^2 + (y(iy)-0.5)^2 + (z(iz)-0.5)^2 )))/0.5);
  % val = exp(-(sum(sqrt( (x(ix)-0)^2 + (y(iy)-0)^2 + (z(iz)-0)^2 )))/0.5);
  % val = sin(4*pi*exp(-(sum(sqrt( (x(ix)-0)^2 + (y(iy)-0)^2 + (z(iz)-0)^2 )))/0.5));
  val = sin(4*pi*exp(-(sum(sqrt( (x(ix)-0)^2 + (y(iy)-0)^2 + (z(iz)-0)^2 )))/0.5)) + sin(4*pi*exp(-(sum(sqrt( (x(ix)-1)^2 + (y(iy)-1)^2 + (z(iz)-1)^2 )))/0.5));

  mat3d(iy,ix,iz) = val;
  mat3d_(iyxz,:)  = [iy,ix,iz,val];
end
% ------------------------------------------------------------------------------
figure;
subplot(1,2,1)
fancy_imagesc(squeeze(mat3d(1,:,:)).')
colormap(rainbow2_cb(1))
colorbar 'off'
axis image;
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
xlabel('x')
ylabel('z')
title('slice of 3d array')
simple_figure()

subplot(1,2,2)
% scatter3(mat3d_(:,2),mat3d_(:,1),mat3d_(:,3),180*ones(ny*nx*nz,1),mat3d_(:,4),'filled')
scatter3(mat3d_(:,2),mat3d_(:,1),mat3d_(:,3),300*abs(mat3d_(:,4)),mat3d_(:,4),'filled')
set(gca,'ZDir','reverse');
axis image;
axis tight;
colormap(rainbow2_cb(1));
colorbar 'off'
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
xlabel('x')
ylabel('y')
zlabel('z')
title('cube of 3d array')
simple_figure()
% ------------------------------------------------------------------------------
figure;

subplot(1,3,1)
scatter3(mat3d_(:,2),mat3d_(:,1),mat3d_(:,3),300*abs(mat3d_(:,4)),mat3d_(:,4),'filled')
set(gca,'ZDir','reverse');
axis image;
axis tight;
colormap(rainbow2_cb(1));
colorbar 'off'
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'Visible','off')
simple_figure()

subplot(1,3,2)
scatter3(mat3d_(:,2),mat3d_(:,1),mat3d_(:,3),300*abs(mat3d_(:,4)),mat3d_(:,4),'filled')
set(gca,'ZDir','reverse');
axis image;
axis tight;
colormap(rainbow2_cb(1));
colorbar 'off'
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'Visible','off')
simple_figure()

subplot(1,3,3)
scatter3(mat3d_(:,2),mat3d_(:,1),mat3d_(:,3),300*abs(mat3d_(:,4)),mat3d_(:,4),'filled')
set(gca,'ZDir','reverse');
axis image;
axis tight;
colormap(rainbow2_cb(1));
colorbar 'off'
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'Visible','off')
simple_figure()
% ------------------------------------------------------------------------------
