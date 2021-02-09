clear
close all
% ------------------------------------------------------------------------------
%  wave
% dt(v) = (1/mu)*grad(u)
% dt(u) = (1/ep)*div(v)
% 
% ------------------------------------------------------------------------------
% --- pde parameters
mu=2;
ep=1;
% --- spatial constraints
X=1;
Z=1;
T=7;
% --- discretization
dx=0.05;
dz=dx;
dt=1e-3;
x=(-1:dx:X).';
z=(-1:dz:Z).';
t=(0:dt:T).';
nx=numel(x);
nz=numel(z);
nt=numel(t);

xx = repmat(x.',nz,1); 
zz = repmat(z,1,nx);
% --- pde parameters
mu=mu*ones(nz,nx);
ep=ep*ones(nz,nx);

mu = ((xx+0.2 )/0.3).^2 + ((zz+0.2)/0.3).^2;
mu = exp(-mu)+1;

ep = ((xx+0.2 )/0.5).^2 + ((zz+0.2)/0.5).^2;
ep = exp(-ep)+0.5;

figure;
subplot(121)
fancy_imagesc(mu,x,z)
colormap(rainbow2(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Magnetic permeability')
simple_figure()

subplot(122)
fancy_imagesc(ep,x,z)
colormap(rainbow2(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Permittivity')
simple_figure()
% ------------------------------------------------------------------------------
% -- space derivative
[Dz_v,Dx_v]=Dx_Dz_v(nz,nx);
[Dz_s,Dx_s]=Dx_Dz_s(nz,nx);
Dz_v = Dz_v/dz;
Dx_v = Dx_v/dx;
Dz_s = Dz_s/dz;
Dx_s = Dx_s/dx;
% -- init
u = zeros(nz,nx);
vx = zeros(nz,nx);
vz = zeros(nz,nx);

u = ((xx-0 )/0.3).^2 + ((zz-0)/0.3).^2;
u = exp(-u);

u_max=max(u(:));
u_min=min(u(:));

u_ = zeros(nz,nx,nt);

clear xx zz
% ------------------------------------------------------------------------------
figure;

subplot(121)
fancy_imagesc(u,x,z)
caxis([u_min u_max])
colormap(rainbow2(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Initial wave')
simple_figure()
% ------------------------------------------------------------------------------
% --- runge kutta
% dt(v) = (1/mu)*grad(u)
% dt(u) = (1/ep)*div(v)
u=u(:);
vx=vx(:);
vz=vz(:);
mu=mu(:);
ep=ep(:);
mu_=1./mu;
ep_=1./ep;
for it=1:(nt-1)
 % - vx
 % runge kuttas
 k1 = mu_.*(Dx_s*u);
 k2 = mu_.*(Dx_s*(u+(dt/2)*k1));
 k3 = mu_.*(Dx_s*(u+(dt/2)*k2));
 k4 = mu_.*(Dx_s*(u+dt*k3));
 vx = vx + (dt/6)*(k1+2*k2+2*k3+k4);
 % - vz
 % runge kuttas
 k1 = mu_.*(Dz_s*u);
 k2 = mu_.*(Dz_s*(u+(dt/2)*k1));
 k3 = mu_.*(Dz_s*(u+(dt/2)*k2));
 k4 = mu_.*(Dz_s*(u+dt*k3));
 vz = vz + (dt/6)*(k1+2*k2+2*k3+k4);
 % - u
 % runge kuttas
 k1 = ep_.*(Dx_v*vx + Dz_v*vz);
 k2 = ep_.*(Dx_v*(vx+(dt/2)*k1) + Dz_v*(vz+(dt/2)*k1));
 k3 = ep_.*(Dx_v*(vx+(dt/2)*k2) + Dz_v*(vz+(dt/2)*k2));
 k4 = ep_.*(Dx_v*(vx+dt*k3) + Dz_v*(vz+dt*k3));
 u = u + (dt/6)*(k1+2*k2+2*k3+k4);
 
 u_(:,:,it) = reshape(u,nz,nx);
end
u=reshape(u,nz,nx);
% ------------------------------------------------------------------------------
% figure;
subplot(122)
fancy_imagesc(u,x,z)
caxis([u_min u_max])
colormap(rainbow2(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Final wave')
simple_figure()
% ------------------------------------------------------------------------------