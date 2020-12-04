clear
close all
% ------------------------------------------------------------------------------
%  Korteweg-de Vries
% dt(u) = - a*u*dx(u) - b*dxxx(u)
% 
% ------------------------------------------------------------------------------
% --- pde parameters
a=1;
b=0.0025;
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
% -- init
u=zeros(nx,nt);
% u(:,1) = exp(-((x-1)/0.8).^2);
% u(:,1) = (sin(x*3));
u(:,1) = cos(x*pi);
% u(:,1) = exp(-((x-1)/0.3).^2)+exp(-((x-0.9)/0.1).^2);
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
% - get 3rd derivative
Dxxx = Dx*(Dx*Dx);
% ------------------------------------------------------------------------------
% --- runge kutta
% dt(u) = - a*u*dx(u) - b*dxxx(u)
for it=1:(nt-1)
% runge kuttas
k1 = -a.*u(:,it).*(Dx*u(:,it)) - b.*Dxxx*u(:,it);
k2 = -a.*(u(:,it)+(dt/2)*k1).*(Dx*(u(:,it)+(dt/2)*k1)) - b.*Dxxx*(u(:,it)+(dt/2)*k1);
k3 = -a.*(u(:,it)+(dt/2)*k2).*(Dx*(u(:,it)+(dt/2)*k2)) - b.*Dxxx*(u(:,it)+(dt/2)*k2);
k4 = -a.*(u(:,it)+dt*k3).*(Dx*(u(:,it)+dt*k3)) - b.*Dxxx*(u(:,it)+dt*k3);

u(:,it+1) = u(:,it) + (dt/6)*(k1+2*k2+2*k3+k4);
end
% ------------------------------------------------------------------------------
figure;
plot(x,u(:,1))
axis tight
xlabel('Length')
ylabel('Units of Korteweg-de Vries')
title('Initial condition')
simple_figure()

figure;
fancy_imagesc(u,t,x)
colormap(rainbow2(1))
axis normal
xlabel('Time')
ylabel('Length')
title('Korteweg-de Vries')
simple_figure()

figure;
fancy_imagesc(u,t,x)
colormap(rainbow2(1))
axis normal
colorbar off
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
