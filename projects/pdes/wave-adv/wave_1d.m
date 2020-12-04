clear
close all
% ------------------------------------------------------------------------------
%  wave advection
% dt(v) = (1/mu)*grad(u)
% dt(u) = (1/ep)*div(v)
% 
% ------------------------------------------------------------------------------
% --- pde parameters
mu=2;
ep=1;
% --- spatial constraints
X=1;
T=1;
% --- discretization
dx=0.005;
dt=1e-5;
x=(-1:dx:X).';
t=(0:dt:T).';
nx=numel(x);
nt=numel(t);
% --- pde parameters
mu=mu*ones(nx,1);
ep=ep*ones(nx,1);

mu=exp(-((x-0.25)/0.3).^2);
mu=mu+1;
ep=exp(-((x-0.25)/0.3).^2);
ep=ep+0.5;
% -- init
u = zeros(nx,nt);
v = zeros(nx,nt);
u(:,1) = exp(-((x-0)/0.3).^2);
% u(:,1) = (sin(x*pi));
% u(:,1) = cos(x*pi);
% u(:,1) = exp(-((x-0)/0.3).^2)+exp(-((x-0.1)/0.1).^2);
% ------------------------------------------------------------------------------
figure;
subplot(311)
plot(x,u(:,1))
axis tight
title('Initial wave')
simple_figure()

subplot(312)
plot(x,mu)
title('Magnetic permeability')
simple_figure()

subplot(313)
plot(x,ep)
xlabel('Length')
title('Permittivity')
simple_figure()
% ------------------------------------------------------------------------------
% -- space derivative
[Dx,~] = Dx_Dz(nx,3);
Dx = Dx(1:nx,1:nx);
Dx = Dx/dx;
% - make periodic boundary conditions!
% Dx(1,1:2)=Dx(2,2:3);
% Dx(1,3)=0;
% Dx(1,nx)=Dx(2,1);
% Dx(nx,(nx-1):nx) = Dx(nx-1,(nx-2):(nx-1));
% Dx(nx,nx-2)=0;
% Dx(nx,1)=Dx(nx-1,nx);
% ------------------------------------------------------------------------------
% --- runge kutta
% dt(v) = (1/mu)*grad(u)
% dt(u) = (1/ep)*div(v)
mu_=1./mu;
ep_=1./ep;
for it=1:(nt-1)
% - v 
% runge kuttas
k1 = mu_.*(Dx*u(:,it));
k2 = mu_.*(Dx*(u(:,it)+(dt/2)*k1));
k3 = mu_.*(Dx*(u(:,it)+(dt/2)*k2));
k4 = mu_.*(Dx*(u(:,it)+dt*k3));

v(:,it+1) = v(:,it) + (dt/6)*(k1+2*k2+2*k3+k4);

% - u
% runge kuttas
k1 = ep_.*(Dx*v(:,it+1));
k2 = ep_.*(Dx*(v(:,it+1)+(dt/2)*k1));
k3 = ep_.*(Dx*(v(:,it+1)+(dt/2)*k2));
k4 = ep_.*(Dx*(v(:,it+1)+dt*k3));

u(:,it+1) = u(:,it) + (dt/6)*(k1+2*k2+2*k3+k4);
end
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(u,t,x)
colormap(rainbow2(1))
axis normal
xlabel('Time')
ylabel('Length')
title('Wave')
simple_figure()

figure;
fancy_imagesc(u,t,x)
colormap(rainbow2(1))
colorbar off
axis normal
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
