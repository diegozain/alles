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
niter=700;
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
rho(fix(nx*0.3):fix(nx*0.5),fix(nz*0.6):fix(nz*0.8)) = 3.5;
rho(fix(nx*0.7):fix(nx*0.8),fix(nz*0.4):fix(nz*0.9)) = 2;
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
ux=ux(:);
uz=uz(:);
% ..............................................................................
% receivers
% ..............................................................................
% these are the fancy ones with 
% 'voltage' type measurements along the z direction
% and placed along the air-ground surface.
% 
% receivers placed every drx in length,
% and at rz_pos and rz_neg in depth.
%
% receivers are placed in the domain with coordinates [rx_ rz_].
% 
% both rx_ and rz_ have the "positive" and "negative" receiver coordinates.
% both rx_ and rz_ consist of two columns: positive and negative.
% each row of rx_ and rz_ make up a receiver "pair".
% ..............................................................................
% drx = 0.5; % m
% rz_pos = 0.1; % m
% rz_neg = 0.6; % m
% % - rx
% rx = 1:drx:(x(nx)-1);
% rx_ = binning(x,rx);
% rx_pos_ = rx_;
% rx_neg_ = rx_;
% rx_ = [rx_pos_ rx_neg_];
% % -rz
% rz_pos_ = binning(z,rz_pos);
% rz_pos_ = repmat(rz_pos_,[numel(rx),1]);
% rz_neg_ = binning(z,rz_neg);
% rz_neg_ = repmat(rz_neg_,[numel(rx),1]);
% rz_= [rz_pos_ rz_neg_];
% ..............................................................................
drx = 0.5; % m
rz_pos = 0.1; % m
rz_neg = 0.6; % m

% - rx
rx = 1:drx:(x(nx)-1);
nr = numel(rx);
rx__ = binning(x,rx);

rx_pos_ = rx__;
rx_neg_ = rx__;
rx_= [rx_pos_ rx_neg_];

% - rz
rz_pos__ = binning(z,rz_pos);
rz_neg__ = binning(z,rz_neg);

rz_pos_ = repmat(rz_pos__,[nr,1]);
rz_neg_ = repmat(rz_neg__,[nr,1]);
rz_= [rz_pos_ rz_neg_];

% - rx + rz
for i_=1:(nr-1)
  rx_pos_ = repmat(rx__(i_),[nr-i_,1]);
  rx_neg_ = rx__((i_+1):nr);
  rx_= [rx_ ; rx_pos_ rx_neg_];
  
  rz_pos_ = repmat(rz_pos__,[nr-i_,1]);
  rz_neg_ = repmat(rz_neg__,[nr-i_,1]);
  rz_= [rz_ ; rz_pos_ rz_neg_];
end
% ..............................................................................
% - rx and rz together
r_ = cat(3,rx_,rz_);
nd = size(r_,1);
r  = zeros(nd,2);
for i_=1:nd
  r(i_,1)=sub2ind([nx,nz],r_(i_,1,1),r_(i_,1,2));
  r(i_,2)=sub2ind([nx,nz],r_(i_,2,1),r_(i_,2,2));
end
% ..............................................................................
figure;
hold on
plot(x(rx_(:,1)),z(rz_(:,1)),'b.','markersize',30)
plot(x(rx_(:,2)),z(rz_(:,2)),'r.','markersize',20)
hold off
legend({'Positive','Negative'})
axis ij
ylim([z(1) z(nz)])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Receivers')
simple_figure()
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
plot(d_x,'r.-','markersize',20)
plot(d_z,'k.-','markersize',20)
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
% - homogeneous background
rho_ = ones(nx,nz);
rho_ = rho_(:);
% - contrast of shapes
[Dx,Dz] = Dx_Dz(nx,nz);
rho_=Dz*rho(:);
rho_(rho_==1.25) = 2;
rho_(rho_==-1.25) = 2;
rho_(rho_==0.5) = 3.5;
rho_(rho_==-0.5) = 3.5;
rho_(rho_==0) = 1;
rho_=reshape(rho_,nx,nz);
ax = 0.9*dx;
az = 0.9*dz;
rho_ = image_gaussian_pad(rho_,az,ax,'LOW_PASS',10,50);
rho_ini = rho_;
rho_=rho_(:);
% ..............................................................................
% 
% inversion
% 
% ..............................................................................
fprintf('\n\n\n     inversion is gonna happen now and will go for %i iterations\n',niter)
% ..............................................................................
% -- optional variables
% store updates
updates=zeros(nz,nx,niter);
% store objective fncs value
Ex=[];
Ez=[];

for ii=1:niter
  % fwd model 
  [ux_,uz_] = g_fwd(Lx,Lz,rho_);

  % data
  d_z_ = M*uz_;
  d_x_ = M*ux_;

  % error (aka residual)
  e_x = d_x_ - d_x;
  e_z = d_z_ - d_z;

  % objective functions value
  Ex = [Ex ; norm(e_x)];
  Ez = [Ez ; norm(e_z)];

  % gradients of objective fnc w/r to density
  g_z = Jzt*e_z;
  g_x = Jxt*e_x;

  % normalize gradients
  g_z = g_z/max(abs(g_z(:)));
  g_x = g_x/max(abs(g_x(:)));

  % kk-filter (smooth gradients)
  ax = 0.9*dx;
  az = 0.9*dz;
  g_z=reshape(g_z,[nx,nz]);
  g_z = image_gaussian_pad(g_z,az,ax,'LOW_PASS',10,50);
  
  ax = 0.9*dx;
  az = 0.9*dz;
  g_x=reshape(g_x,[nx,nz]);
  g_x = image_gaussian_pad(g_x,az,ax,'LOW_PASS',10,50);
  
  g_z=g_z(:)/max(abs(g_z(:)));
  g_x=g_x(:)/max(abs(g_x(:)));
  
  g_x=g_x+g_z;
  g_z=g_x;

  % step size
  k=5e-2; % 1e-1
  drho = k*g_z;
  [step_x,step_z] = g_pica(d_x,d_z,e_x,e_z,M,Lx,Lz,rho_,drho,k);
  % update (noun)
  % update = step_z*g_z;
  % update = step_x*g_x;
  update = 0.5*(step_x*g_x + step_z*g_z);
  % update (verb)
  rho_ = rho_ .* exp( -rho_.*update );
  
  % store updates
  updates(:,:,ii) = reshape(update,nx,nz).';
  
  if mod(ii,fix(0.3*niter))==0
    fprintf('\n just finished iteration #%i',ii);
  end
end
fprintf('\n\n')
rho_=reshape(rho_,[nx,nz]);
% ..............................................................................
mini = min([min(rho(:)) min(rho_(:))]);
maxi = max([max(rho(:)) max(rho_(:))]);
% ..............................................................................
figure;
hold on;
plot((1:niter),Ex,'b.-','markersize',15);
plot((1:niter),Ez,'r.-','markersize',15);
hold off;
legend({'x','z'})
xlabel('Iteration #')
title('Objective function values')
simple_figure()
% ..............................................................................
figure('Renderer', 'painters', 'Position', [10 10 500 400]);
subplot(3,1,1)
hold on;
fancy_imagesc(reshape(rho,[nx,nz]).',x,z)
colorbar('off')
plot(x(rx_(:,1)),z(rz_(:,1)),'k.','markersize',25)
plot(x(rx_(:,2)),z(rz_(:,2)),'w.','markersize',15);axis ij
% plot(rx,zeros(size(rx)),'k.','markersize',15);axis ij
hold off;
caxis([mini maxi])
axis image
% ylabel('Depth (m)')
% xlabel('Length (m)')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('True density')
simple_figure()

subplot(3,1,2)
fancy_imagesc(rho_ini.',x,z)
colorbar('off')
caxis([mini maxi])
axis image
% xlabel('Length (m)')
% ylabel('Depth (m)')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Initial density')
simple_figure()

subplot(3,1,3)
fancy_imagesc(rho_.',x,z)
colorbar('off')
caxis([mini maxi])
axis image
% xlabel('Length (m)')
% ylabel('Depth (m)')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Recovered density')
simple_figure()
% ..............................................................................
