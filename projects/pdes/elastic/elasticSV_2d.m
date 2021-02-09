clear
close all
% ------------------------------------------------------------------------------
%  2d elastic SV wave
% 
% rho*dt(v) = grad * T + F
% s*dt(u) = grad.' * v - tau*T
%
% where,
% 
% v = [vx ; vz]       particle velocity
% T = [sxx; sxz; szz] stress
% 
% grad = [dx 0  dz]
%        [0  dz dx]
% 
% inv(s) = [lam+2*mu lam        0]
%           [0         0        mu]
%           [lam       lam+2*mu  0]
% 
% ------------------------------------------------------------------------------
% vp = sqrt((lam + 2*mu)./rho);
% vs = sqrt(mu./rho);
% 
% lam = (vp.^2-2*vs.^2).*rho;
% mu  = rho.*vs^2;
% ------------------------------------------------------------------------------
% --- pde parameters
lam_= [1e6 ; 4e7]; % [1; 2];
mu_ = [1e7 ; 3e8]; % [2; 4];
rho_= [1.7e3 ; 2e3]; % [3; 6];
% --- spatial constraints
X=2;
Z=2;
T=15;
% --- source parameters
to = 4;
fo = 1;
wo = 2*pi*fo;
% --- stability
vp = sqrt((lam_ + 2*mu_)./rho_);
vs = sqrt(mu_./rho_);

vs_min = min(vs);
vp_min = min(vp);
vs_max = max(vs);
vp_max = max(vp);

if vs==0
  vel_min = vp_min;
else
  vel_min = min([vs_min,vp_min]);
end
vel_max = max([vs_max,vp_max]);

fmax = 2.2*fo;
% - space
l_min = vel_min / fmax;
no_p_wa = 10; % 10 50
dx = l_min/no_p_wa;
dz = dx;
% - time
courant_factor = 0.9;
dt = 1/(vel_max * sqrt((1/dx^2)+(1/dz^2)));
dt = courant_factor * dt;
% --- discretization
x=(-X:dx:X).';
z=(-X:dz:Z).';
t=(0:dt:T).';
nx=numel(x);
nz=numel(z);
nt=numel(t);
% --- pde parameters
lam=lam_(1)*ones(nz,nx);
mu =mu_(1)*ones(nz,nx);
rho=rho_(1)*ones(nz,nx);

lam(fix(nz*(1/3)):fix(nz*(2/3)),fix(nx*(1/3)):fix(nx*(2/3))) = lam_(2);
mu(fix(nz*(1/3)):fix(nz*(2/3)),fix(nx*(1/3)):fix(nx*(2/3))) = mu_(2);
rho(fix(nz*(1/3)):fix(nz*(2/3)),fix(nx*(1/3)):fix(nx*(2/3))) = rho_(2);

lam2mu = lam+2*mu;

figure;
subplot(131)
fancy_imagesc(lam,x,z)
colormap(rainbow2(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Lame #1')
simple_figure()

subplot(132)
fancy_imagesc(mu,x,z)
colormap(rainbow2(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Lame #2')
simple_figure()

subplot(133)
fancy_imagesc(rho,x,z)
colormap(rainbow2(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Density')
simple_figure()
% ------------------------------------------------------------------------------
% -- source function
src_xz = [(x(end)+x(1))/2 , (z(end)+z(1))/2];
src_ix = binning(x,src_xz(1));
src_iz = binning(z,src_xz(2));
isrc = sub2ind([nz,nx],src_iz,src_ix);

fx=zeros(nt,1);
fz=( 1-0.5*(wo^2)*(t-to).^2 ) .* exp( -0.25*(wo^2)*(t-to).^2 );

figure;
plot(t,fz)
xlabel('Time')
ylabel('Amplitude')
title('Source f_z')
simple_figure()
% ------------------------------------------------------------------------------
% -- space derivative
[Dz_v,Dx_v]=Dx_Dz_v(nz,nx);
[Dz_s,Dx_s]=Dx_Dz_s(nz,nx);
% -- init
vx = zeros(nz,nx);
vz = zeros(nz,nx);

sxx=zeros(nz,nx);
szz=zeros(nz,nx);
sxz=zeros(nz,nx);

vz_=zeros(nz,nx,nt);

rho_=1./rho;
% ------------------------------------------------------------------------------
% % --- wave solver
% vx=vx(:);
% vz=vz(:);
% 
% sxx=sxx(:);
% szz=szz(:);
% sxz=sxz(:);
% 
% mu=mu(:);
% lam=lam(:);
% lam2mu=lam2mu(:);
% rho_=rho_(:);
% 
% tic;
% for it=1:nt
%  % - velocity
% 
%  dxsxx_dzsxz = (1/dx)*Dx_v*sxx + (1/dz)*Dz_v*sxz;
%  dxsxz_dzszz = (1/dx)*Dx_v*sxz + (1/dz)*Dz_v*szz;
% 
%  % - source terms
%  dxsxx_dzsxz(isrc) = dxsxx_dzsxz(isrc) + fx(it);
%  dxsxz_dzszz(isrc) = dxsxz_dzszz(isrc) + fz(it);
% 
%  vx = vx + dt*rho_.*dxsxx_dzsxz;
%  vz = vz + dt*rho_.*dxsxz_dzszz;
% 
%  % - absorbing boundary conditions
% 
%  % - store wavefield
%  vz_(:,:,it) = reshape(vz,nz,nx);
% 
%  % - stress
%  dxvx = (1/dx)*Dx_s*vx;
%  dxvz = (1/dx)*Dx_s*vz;
%  dzvx = (1/dz)*Dz_s*vx;
%  dzvz = (1/dz)*Dz_s*vz;
% 
%  sxx=sxx + dt*( (lam2mu).*dxvx + lam.*dzvz );
%  szz=szz + dt*( (lam2mu).*dzvz + lam.*dxvx );
%  sxz=sxz + dt*mu.*( dzvx + dxvz );
% end
% toc;
% 
% vx=reshape(vx,nz,nx);
% vz=reshape(vz,nz,nx);
% 
% sxx=reshape(sxx,nz,nx);
% szz=reshape(szz,nz,nx);
% sxz=reshape(sxz,nz,nx);
% ------------------------------------------------------------------------------
tic;
for it=1:nt
  % -- velocity updates
  dxsxx_dzsxz=zeros(nz,nx);
  dxsxz_dzszz=zeros(nz,nx);

  % 2nd order
  ix=2:nx;
  % dx(sxx)
  dxsxx_dzsxz(:,ix)=(sxx(:,ix)-sxx(:,ix-1))/dx;
  % dx(sxz)
  dxsxz_dzszz(:,ix)=(sxz(:,ix)-sxz(:,ix-1))/dx;
  iz=2:nz;
  % dz(sxz)
  dxsxx_dzsxz(iz,:)=dxsxx_dzsxz(iz,:) + (sxz(iz,:)-sxz(iz-1,:))/dz;
  % dz(szz)
  dxsxz_dzszz(iz,:)=dxsxz_dzszz(iz,:) + (szz(iz,:)-szz(iz-1,:))/dz;

  % % 4th order
  % for ix=3:nx-1
  %  % dx(sxx)
  %  dxsxx_dzsxz(:,ix)=(9/8)*(sxx(:,ix)-sxx(:,ix-1))/dx - (1/24)*(sxx(:,ix+1)-sxx(:,ix-2))/dx;
  %  % dx(sxz)
  %  dxsxz_dzszz(:,ix)=(9/8)*(sxz(:,ix)-sxz(:,ix-1))/dx - (1/24)*(sxz(:,ix+1)-sxz(:,ix-2))/dx;
  % end
  % for iz=3:nz-1
  %  % dz(sxz)
  %  dxsxx_dzsxz(iz,:)=dxsxx_dzsxz(iz,:) + (9/8)*(sxz(iz,:)-sxz(iz-1,:))/dz - (1/24)*(sxz(iz+1,:)-sxz(iz-2,:))/dz;
  %  % dz(szz)
  %  dxsxz_dzszz(iz,:)=dxsxz_dzszz(iz,:) + (9/8)*(szz(iz,:)-szz(iz-1,:))/dz - (1/24)*(szz(iz+1,:)-szz(iz-2,:))/dz;
  % end

  % % different grid?
  % dxsxx=zeros(nz,nx);
  % dxsxz=zeros(nz,nx);
  % dzsxz=zeros(nz,nx);
  % dzszz=zeros(nz,nx);
  % ix=3:nx-1;
  % dxsxx(:,ix)=(9/8)*(sxx(:,ix)-sxx(:,ix-1))/dx - (1/24)*(sxx(:,ix+1)-sxx(:,ix-2))/dx;
  % iz=2:nz-2;
  % dzsxz(iz,:)=(9/8)*(sxz(iz+1,:)-sxz(iz,:))/dz - (1/24)*(sxz(iz+2,:)-sxz(iz-1,:))/dz;
  % dxsxx_dzsxz =  dxsxx+dzsxz;
  % iz=3:nz-1;
  % dzszz(iz,:) = (9/8)*(szz(iz,:)-szz(iz-1,:))/dz - (1/24)*(szz(iz+1,:)-szz(iz-2,:))/dz;
  % ix=2:nx-2;
  % dxsxz(:,ix) = (9/8)*(sxz(:,ix+1)-sxz(:,ix))/dx - (1/24)*(sxz(:,ix+2)-sxz(:,ix-1))/dx;
  % dxsxz_dzszz = dxsxz+dzszz;

  % - velocity source 
  dxsxx_dzsxz(isrc) = dxsxx_dzsxz(isrc) + fx(it);
  dxsxz_dzszz(isrc) = dxsxz_dzszz(isrc) + fz(it);

  vx = vx + dt*rho_.*dxsxx_dzsxz;
  vz = vz + dt*rho_.*dxsxz_dzszz;
  
  % - absorbing boundary conditions
  

  % - store wavefield
  vz_(:,:,it) = reshape(vz,nz,nx);

  % -- stress updates
  dxvx=zeros(nz,nx);
  dzvx=zeros(nz,nx);
  dxvz=zeros(nz,nx);
  dzvz=zeros(nz,nx);

  % 2nd order
  ix=1:nx-1;
  % dx(vx)
  dxvx(:,ix)=(vx(:,ix+1)-vx(:,ix))/dx;
  % dx(vz)
  dxvz(:,ix)=(vz(:,ix+1)-vz(:,ix))/dx;
  iz=1:nz-1;
  % dz(vx)
  dzvx(iz,:)=(vx(iz+1,:)-vx(iz,:))/dz;
  % dz(vz)
  dzvz(iz,:)=(vz(iz+1,:)-vz(iz,:))/dz;

  % % 4th order
  % for ix=2:nx-2
  %  % dx(vx)
  %  dxvx(:,ix)=9*(vx(:,ix+1)-vx(:,ix))/(8*dx)-(vx(:,ix+2)-vx(:,ix-1))/(24*dx);
  %  % dx(vz)
  %  dxvz(:,ix)=9*(vz(:,ix+1)-vz(:,ix))/(8*dx)-(vz(:,ix+2)-vz(:,ix-1))/(24*dx);
  % end
  % for iz=2:nz-2
  %  % dz(vx)
  %  dzvx(iz,:)=9*(vx(iz+1,:)-vx(iz,:))/(8*dz)-(vx(iz+2,:)-vx(iz-1,:))/(24*dz);
  %  % dz(vz)
  %  dzvz(iz,:)=9*(vz(iz+1,:)-vz(iz,:))/(8*dz)-(vz(iz+2,:)-vz(iz-1,:))/(24*dz);
  % end

  % % different grid?
  % ix=2:nx-2;
  % dxvx(:,ix)=9*(vx(:,ix+1)-vx(:,ix))/(8*dx)-(vx(:,ix+2)-vx(:,ix-1))/(24*dx);
  % iz=2:nz-2;
  % dzvz(iz,:)=9*(vz(iz+1,:)-vz(iz,:))/(8*dz)-(vz(iz+2,:)-vz(iz-1,:))/(24*dz);
  % iz=3:nz-1;
  % dzvx(iz,:)=9*(vx(iz,:)-vx(iz-1,:))/(8*dz)-(vx(iz+1,:)-vx(iz-2,:))/(24*dz);
  % ix=3:nx-1;
  % dxvz(:,ix)=9*(vz(:,ix)-vz(:,ix-1))/(8*dx)-(vz(:,ix+1)-vz(:,ix-2))/(24*dx);

  sxx=sxx + dt*( (lam2mu).*dxvx + lam.*dzvz );
  szz=szz + dt*( (lam2mu).*dzvz + lam.*dxvx );
  sxz=sxz + dt*mu.*( dzvx + dxvz );
end
toc;
% ------------------------------------------------------------------------------
vz_min = min(vz_(:));
vz_max = max(vz_(:));

figure;
fancy_imagesc(vz,x,z)
colormap(rainbow2(1))
caxis([vz_min vz_max])
hold on
plot(src_xz(1),src_xz(2),'c.','markersize',50)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Final v_z')
simple_figure()
% ------------------------------------------------------------------------------