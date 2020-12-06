close all
clear
clc
% ..............................................................................
addpath('src/')
% ..............................................................................
% diego domenzain
% fall 2017
% Boise State University
% ..............................................................................
max_iter=200;
% ..............................................................................
fprintf('\n     given gravimetry data on the x and z direction\n     only at surface receivers,\n     recover density.')

fprintf('\n\n         receivers are shown in the last plot as solid dots\n')
% ..............................................................................
% gravity
dx=0.1;
dz=0.1;
x = (0:dx:12).';
z = (0:dz:3).';
nx = numel(x);
nz = numel(z);
% build geometry dependent linear operators
fprintf('\n\n         ~-~-~ building linear operators ~-~-~\n')
[Lx,Lz] = g_L(x,z);
% ..............................................................................
% true density
% ..............................................................................
% bg
rho = ones(nx,nz);
% box
rho(fix(nx*0.4):fix(nx*0.6),fix(nz*0.4):fix(nz*0.6)) = 3.5;
% % smooth blob
% x_ = linspace(-1,1,nx);
% z_ = linspace(-1,1,nz);
% [Y,X] = meshgrid(z_,x_);
% sigma = 0.2;
% rho = 1 + 8 * exp(-((X-0.5).^2+(Y).^2)/(2*sigma^2));
figure;
fancy_imagesc(rho.',x,z)
axis image
xlabel('Length (m)')
ylabel('Depth (m)')
title('True density')
simple_figure()
% ..............................................................................
% make column
rho = rho(:);
% ..............................................................................
% 
% fwd gravity
% 
% ..............................................................................
% gravity x,z
fprintf('\n\n         *-*-* computing synthetic data *-*-*\n')
[ux,uz] = g_fwd(Lx,Lz,rho);
% ..............................................................................
% see grav field
uz=reshape(uz,[nx,nz]);
ux=reshape(ux,[nx,nz]);
% ..............................................................................
figure;
subplot(2,1,1)
fancy_imagesc(ux.',x,z)
ylabel('Depth (m)')
xlabel('Length (m)')
title('Component x of gravity' )
simple_figure()
subplot(2,1,2)
fancy_imagesc(uz.',x,z)
xlabel('Length (m)')
ylabel('Depth (m)')
title('Component z of gravity' )
simple_figure()
% ..............................................................................
% make column
ux=ux(:);uz=uz(:);
% ..............................................................................
% receivers
% ..............................................................................
% these are the fancy ones with 
% 'voltage' type measurements along the z direction
% and placed along the air-ground surface (between 0.5m and 1m depth).
rx = 1:0.5:(x(end)-1);
rx_ = binning(x,rx);
rx_ = [rx_ rx_];
rz_pos = 0.5;
rz_pos_ = binning(z,rz_pos);
rz_pos_ = repmat(rz_pos_,[numel(rx),1]);
rz_neg = 1;
rz_neg_ = binning(z,rz_neg);
rz_neg_ = repmat(rz_neg_,[numel(rx),1]);
rz_ = [rz_pos_ rz_neg_];
r_=cat(3,rx_,rz_);
nd=size(r_,1);
r=zeros(nd,2);
for i_=1:nd
r(i_,1)=sub2ind([nx,nz],r_(i_,1,1),r_(i_,1,2));
r(i_,2)=sub2ind([nx,nz],r_(i_,2,1),r_(i_,2,2));
end
% ..............................................................................
% measuring operator
M = g_M(x,z,r);
% ..............................................................................
% data
d_z = M*uz;
d_x = M*ux;
% ..............................................................................
% see data
figure;
hold on
plot(rx,d_x,'r.-','markersize',20)
plot(rx,d_z,'k.-','markersize',20)
hold off
axis tight
legend({'\Delta u_x','\Delta u_z'})
xlabel('Length (m)')
title('Data' )
simple_figure()
% ..............................................................................
% jacobians of data w/r to density
Jzt=(M*Lz).';
Jxt=(M*Lx).';
% ..............................................................................
% initial guess
rho_ = ones(nx,nz);
rho_ = rho_(:);
% ..............................................................................
% 
% inversion
% 
% ..............................................................................
fprintf('\n\n\n     inversion is gonna happen now and will go for %i iterations\n',max_iter)
% ..............................................................................
for ii=1:max_iter
  % fwd model 
  [ux_,uz_] = g_fwd(Lx,Lz,rho_);

  % data
  d_z_ = M*uz_;
  d_x_ = M*ux_;

  % error (aka residual)
  e_x = d_x_ - d_x;
  e_z = d_z_ - d_z;

  % gradients of objective fnc w/r to density
  g_z = Jzt*e_z;
  g_x = Jxt*e_x;

  % normalize gradients
  g_z = g_z/max(abs(g_z(:)));
  g_x = g_x/max(abs(g_x(:)));

  % kk-filter (smooth gradients)
  ax = 0.8*dx; az = 0.8*dz;
  g_z=reshape(g_z,[nx,nz]);
  [g_z, g_z_kk, ~] = image_gaussian(g_z,az,ax,'LOW_PASS');
  g_z=g_z(:)/max(abs(g_z(:)));
  
  ax = 0.5*dx; az = 0.5*dz;
  g_x=reshape(g_x,[nx,nz]);
  [g_x, g_x_kk, ~] = image_gaussian(g_x,az,ax,'LOW_PASS');
  g_x=g_x(:)/max(abs(g_x(:)));

  % step size
  k=0.1;
  drho = k*g_z;
  [step_x,step_z] = g_pica(d_x,d_z,e_x,e_z,M,Lx,Lz,rho_,drho,k);
  % update
  rho_ = rho_ .* exp( -rho_.*(step_z*g_z) );
  
  if mod(ii,fix(0.3*max_iter))==0
    fprintf('\n just finished iteration #%i',ii);
  end
end
fprintf('\n\n')
rho_=reshape(rho_,[nx,nz]);
% ..............................................................................
mini = min([min(rho(:)) min(rho_(:))]);
maxi = max([max(rho(:)) max(rho_(:))]);
% ..............................................................................
figure('Renderer', 'painters', 'Position', [10 10 500 400]);
subplot(2,1,1)
hold on;
fancy_imagesc(reshape(rho,[nx,nz]).',x,z)
colorbar('off')
plot(rx,zeros(size(rx)),'k.','markersize',15);axis ij
hold off;
caxis([mini maxi])
axis image
ylabel('Depth (m)')
xlabel('Length (m)')
title('True density')
simple_figure()

subplot(2,1,2)
fancy_imagesc(rho_.',x,z)
colorbar('off')
caxis([mini maxi])
axis image
xlabel('Length (m)')
ylabel('Depth (m)')
title('Recovered density')
simple_figure()
% ..............................................................................
