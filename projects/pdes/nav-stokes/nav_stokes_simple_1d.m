clear
close all
% ------------------------------------------------------------------------------
% 
% dt(h) = - sqrt(So)/n * dx(h^(5/3)) + q
%
% this is a simplified navier stokes equation that models the 
% height (h) of a fluid over time on a plane that has a slope controlled by So.
%
% the navier stokes simplification is taken from god knows where.
% 
% ------------------------------------------------------------------------------
% --- pde parameters
n=0.025;
So=0.02;
r=2.3e-5;
b=sqrt(So)/n;
% --- spatial constraints
X=1;
T=1.5;
% --- discretization
dx=1e-3; % 4.5
dt=1e-4;  % 0.5
x=0:dx:X;
t=0:dt:T;
nx=numel(x);
nt=numel(t);
% --- source term
q=r*t;
d=binning(t,t(nt)*0.3);
q(d:nt)=0;
% q=exp(-((t-t(d))/0.05).^2);
% ------------------------------------------------------------------------------
% --- runge kutta v2
% dt(h) = sqrt(So)/n * dx(h^(5/3)) + q
% dt(h) = b * dx(h^(5/3)) + q
% dx(h^(5/3)) = (5/3)*h^(2/3) * dx(h)
% dt(h) = b * ((5/3)*h^(2/3) * dx(h)) + q
% -- space derivative
[Dx,~] = Dx_Dz(nx,3);
Dx = Dx(1:nx,1:nx);
Dx = Dx/dx;
% -- init
h = zeros(nx,nt);
for it=1:(nt-1)
% runge kuttas
% k1 = b.*(Dx*(h(:,it).^(5/3)));
% k2 = b.*( Dx*( ( h(:,it) + (dt/2)*k1 ).^(5/3) ) );
% k3 = b.*( Dx*( ( h(:,it) + (dt/2)*k2 ).^(5/3) ) );
% k4 = b.*( Dx*( ( h(:,it) + dt*k3 ).^(5/3) ) );

k1 = b.*( ((5/3)*h(:,it).^(2/3)) .* (Dx*h(:,it)) );
k2 = b.*( ((5/3)*(h(:,it)+ (dt/2)*k1).^(2/3)) .* (Dx*(h(:,it)+ (dt/2)*k1)) );
k3 = b.*( ((5/3)*(h(:,it)+ (dt/2)*k2).^(2/3)) .* (Dx*(h(:,it)+ (dt/2)*k2)) );
k4 = b.*( ((5/3)*(h(:,it)+ dt*k3).^(2/3)) .* (Dx*(h(:,it)+ dt*k3)) );

h(:,it+1) = h(:,it) + (dt/6)*(k1+2*k2+2*k3+k4) + q(it);

h(nx,it+1) = 0;
end
h=h.';
% ------------------------------------------------------------------------------
figure;
subplot(211)
plot(t,q)
axis tight
xlabel('Time (s)')
ylabel('Source ( ? )')
title('Source')
simple_figure()
subplot(212)
plot(t,h(:,binning(x,0.5*x(nx))))
axis tight
xlabel('Time (s)')
ylabel('h (m)')
title('Midpoint')
simple_figure()
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(h,x,t)
colormap(parula)
axis normal
ylabel('Time (s)')
xlabel('Length (m)')
title('Height of fluid (m)')
simple_figure()
% ------------------------------------------------------------------------------
fprintf('\n\n     animation time!!\n')
fprintf('\n\n       press any key multiple times to see\n       the height advance through time.\n')
figure;
for i_=1:100:nt
 plot(x,h(i_,:))
 xlabel('Length (m)')
 ylabel('h (m)')
 ylim([min(h(:)) max(h(:))])
 simple_figure()
 pause
end
% ------------------------------------------------------------------------------
