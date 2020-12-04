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
k=2;
a=3;
% --- spatial constraints
X=1;
T=0.1;
% --- discretization
dx=0.005;
dt=1e-6;
x=(-1:dx:X).';
t=(0:dt:T).';
nx=numel(x);
nt=numel(t);
% --- pde parameters
k=k*ones(nx,1);
a=a*ones(nx,1);

k=exp(-((x-0.25)/0.3).^2);
k=k+1;
a=exp(-((x-0.25)/0.3).^2);
a=a+3;
% -- init
u = zeros(nx,nt);
% u(:,1) = exp(-((x-0)/0.3).^2);
% u(:,1) = (sin(x*pi));
% u(:,1) = cos(x*pi);
u(:,1) = exp(-((x+0.4)/0.1).^2)+exp(-((x-0.4)/0.1).^2);
% ------------------------------------------------------------------------------
figure;
subplot(311)
plot(x,u(:,1))
axis tight
title('Initial heat')
simple_figure()

subplot(312)
plot(x,k)
title('Kappa')
simple_figure()

subplot(313)
plot(x,a)
xlabel('Length')
title('Thermal conductivity')
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
% --
% Dx(1,1:3) = 0;
% Dx(nx,(nx-2):nx) = 0;
% ------------------------------------------------------------------------------
% --- runge kutta
% dt(u) = k div a grad u
a_x=Dx*a;
Dxx=Dx*Dx;
for it=1:(nt-1)
% runge kuttas
k1 = k.*( a_x.*(Dx*u(:,it)) + a.*(Dxx*u(:,it)) );
k2 = k.*( a_x.*(Dx*(u(:,it)+(dt/2)*k1)) + a.*(Dxx*(u(:,it)+(dt/2)*k1)) );
k3 = k.*( a_x.*(Dx*(u(:,it)+(dt/2)*k2)) + a.*(Dxx*(u(:,it)+(dt/2)*k2)) );
k4 = k.*( a_x.*(Dx*(u(:,it)+dt*k3)) + a.*(Dxx*(u(:,it)+dt*k3)) );

u(:,it+1) = u(:,it) + (dt/6)*(k1+2*k2+2*k3+k4);
end
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(u,t,x)
colormap(rainbow2(1))
axis normal
xlabel('Time')
ylabel('Length')
title('Final heat')
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
