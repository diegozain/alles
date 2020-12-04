clear
close all
% ------------------------------------------------------------------------------
%  heat
% 
% dt(u) = k div a grad u
%
%  The coefficient k(x) is the inverse of specific heat of the substance at x 
%  times density of the substance at x: k= 1 / ( rho * c_p ).
% ------------------------------------------------------------------------------
% --- pde parameters
k=1;
a=4;
% --- spatial constraints
X=1;
Z=1;
T=5;
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
k=k*ones(nz,nx);
a=a*ones(nz,nx);

k = ((xx+0.2 )/0.3).^2 + ((zz+0.2)/0.3).^2;
k = 2*exp(-k)+1;

a = ((xx+0.2 )/0.5).^2 + ((zz+0.2)/0.5).^2;
a = 5*exp(-a)+1;

figure;
subplot(121)
fancy_imagesc(k,x,z)
colormap(rainbow2(1))
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
colorbar off
title('Kappa')
simple_figure()

subplot(122)
fancy_imagesc(a,x,z)
colormap(rainbow2(1))
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
colorbar off
title('Thermal conductivity')
simple_figure()
% ------------------------------------------------------------------------------
% -- space derivative
[Dz,Dx] = Dx_Dz(nz,nx);
% -- init
u = zeros(nz,nx);

u = ((xx-0 )/0.3).^2 + ((zz-0)/0.3).^2;
u = exp(-u);

u_max=max(u(:));
u_min=min(u(:));

clear xx zz
% ------------------------------------------------------------------------------
figure;

subplot(121)
fancy_imagesc(u,x,z)
caxis([u_min u_max])
colormap(rainbow2(1))
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
colorbar off
title('Initial heat')
simple_figure()
% ------------------------------------------------------------------------------
% --- runge kutta
% dt(u) = k div a grad u
u=u(:);
a=a(:);
k=k(:);
a_x=Dx*a;
a_z=Dz*a;
Dxx=Dx*Dx;
Dzz=Dz*Dz;
for it=1:(nt-1)
% runge kuttas
k1 = k.*( a_x.*(Dx*u) + a.*(Dxx*u) + a_z.*(Dz*u) + a.*(Dzz*u));
k2 = k.*( a_x.*(Dx*(u+(dt/2)*k1)) + a.*(Dxx*(u+(dt/2)*k1)) + a_z.*(Dz*(u+(dt/2)*k1)) + a.*(Dzz*(u+(dt/2)*k1)));
k3 = k.*( a_x.*(Dx*(u+(dt/2)*k2)) + a.*(Dxx*(u+(dt/2)*k2)) + a_z.*(Dz*(u+(dt/2)*k2)) + a.*(Dzz*(u+(dt/2)*k2)));
k4 = k.*( a_x.*(Dx*(u+dt*k3)) + a.*(Dxx*(u+dt*k3)) + a_z.*(Dz*(u+dt*k3)) + a.*(Dzz*(u+dt*k3)));

u = u + (dt/6)*(k1+2*k2+2*k3+k4);
end
u=reshape(u,nz,nx);
% ------------------------------------------------------------------------------
% figure;
subplot(122)
fancy_imagesc(u,x,z)
caxis([u_min u_max])
colormap(rainbow2(1))
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
colorbar off
title('Final heat')
simple_figure()
% ------------------------------------------------------------------------------