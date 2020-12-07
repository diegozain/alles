clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
%          eikonal solver
% ------------------------------------------------------------------------------
nx = 100; 
nz = 50;
% ------------------------------------------------------------------------------
% adjacency matrix
A = adj_grid(nx,nz);
% ------------------------------------------------------------------------------
% velocity on nodes
c = ones(nx,nz);

% % a rectangle in the middle
% c(fix(nx/3):fix(2*nx/3),fix(nz/3):fix(nz/2)) = 3.5;
% c(fix(nx/3):fix(2*nx/3),fix(nz/2):fix(2*nz/3)) = 3.5;

% a gaussian blob in the middle
x = linspace(-1,1,nx);
z = linspace(-1,1,nz);
[Y,X] = meshgrid(z,x);
sigma = 0.5;
c = 1 + 8 * exp(-(X.^2+Y.^2)/(2*sigma^2));
% ------------------------------------------------------------------------------
% choose source node
% source = nx*z + x
s = [nx*5 + ceil(nx*0.1)];
% ------------------------------------------------------------------------------
h  = 0.1;
tic;
[u,P] = fast_marcher(A,c,s,h);
toc;
% ------------------------------------------------------------------------------
u = reshape(u,[nx,nz]);
% ------------------------------------------------------------------------------
% geodesic paths
% ------------------------------------------------------------------------------
% choose receiver
% receiver = nx*z + x
r = nx*5 + ceil(nx*0.9);
% gradient of u
grad_u = grad_normal(u);
% ray path
p = ray_path(r,s,grad_u,h,nx);
% ------------------------------------------------------------------------------
% bring to matrix notation
r = p2ij(r,nx);
s = p2ij(s,nx);
% ------------------------------------------------------------------------------
%       see
% ------------------------------------------------------------------------------
% % this one just looks cool, but it's not anything physical
% displ = @(v)cos(2*pi*5*v/max(v(:)));
% 
% figure;
% fancy_imagesc(displ(u.'),1:nx,1:nz);
% colorbar('off')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% simple_figure();
% ------------------------------------------------------------------------------
figure('Renderer', 'painters', 'Position', [10 10 700 200]);

subplot(121)
fancy_imagesc(u.',1:nx,1:nz);
colorbar('off')
hold on
plot(p(:,1),p(:,2),'k-')
plot(r(:,1),r(:,2),'.','MarkerSize',40,'color',[0.1412,0.8471,0.1882])
plot(s(:,1),s(:,2),'rp','MarkerSize',20,'MarkerFaceColor','r')
hold off
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Eikonal wavefield')
simple_figure();

subplot(122)
fancy_imagesc(c.',1:nx,1:nz);
colorbar('off')
hold on
plot(p(:,1),p(:,2),'k-')
plot(r(:,1),r(:,2),'m.','MarkerSize',40,'color',[0.1412,0.8471,0.1882])
plot(s(:,1),s(:,2),'rp','MarkerSize',20,'MarkerFaceColor','r')
hold off
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Velocity')
simple_figure();
% ------------------------------------------------------------------------------
