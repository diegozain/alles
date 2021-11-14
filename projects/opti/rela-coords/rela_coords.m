clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% given distances among 4 points,
% return their coordinates.
%
% distances are given along the perimeter & diagonals:
%
%       1
%    2.----.1
%     |\5 /|
%    2| /\ |4
%     |/3 \|
%    3.-----.4
%        6
% ------------------------------------------------------------------------------
% observed data
r = [4.5; 3.0; 6.5; 4.8; 3.6; 3.85];
r = r.^2;
% initial guess
% x = [6; 3; -1; 0; 4; 4];
x = [4.8; 3; -1; 0; 4; 4];
% x = 1e1*randn(6,1);
% ------------------------------------------------------------------------------
figure(1);
hold on;
plot(x(1:3),x(4:6),'r.','markersize',25);
plot(0,0,'r.','markersize',25)
hold off;
xlabel('Length (m)')
ylabel('Width (m)')
title('Initial & final')
simple_figure()
% ------------------------------------------------------------------------------
gradd_ = zeros(6,1);
% ------------------------------------------------------------------------------
niter  = 100;
OBJ = 'rms_log';
objfnc_= Inf;
objfnc = zeros(niter,1);
tol_error = -11;%1e-5;
% ------------------------------------------------------------------------------
iter   = 1;
while (objfnc_>tol_error & iter<=niter)
  % % 🎨
  % figure(1001);
  % hold on;
  % plot(x(1:3),x(4:6),'.','markersize',20)
  % plot(0,0,'.','markersize',20)
  % hold off;
  % axis image;
  % xlabel('Length (m)')
  % ylabel('Width (m)')
  % simple_figure()

  % forward model 👉
  dat_ = fwd_dista(x);

  % obj function 👇
  [objfnc_,err_]= obj_dista(dat_, r, OBJ);
  objfnc(iter)  = objfnc_;

  % bacwkard model 👈
  Jac  = jac_dista(x);
  gra_ = Jac * err_;

  % % this should only be done when the first point is also fixed in
  % % the initial condition.
  % gra_(1)=0;gra_(4)=0;

  % step size 🚶
  step_= step_dista(x,r,gra_,iter,OBJ);
  gra_ = step_*gra_;

  % update
  x    = x - gra_;
  iter = iter+1;
end
% ------------------------------------------------------------------------------
%
%                             🎨🎨🎨🎨
%
% ------------------------------------------------------------------------------
figure(1);
hold on;
plot(x(1:3),x(4:6),'k.','markersize',25)
plot(0,0,'k.','markersize',25)
hold off;
axis image;
xlabel('Length (m)')
ylabel('Width (m)')
title('Final')
simple_figure()
% ------------------------------------------------------------------------------
figure;
hold on;
plot(x(1:3),x(4:6),'k.','markersize',30)
plot(0,0,'k.','markersize',30)
hold off;
grid on;
axis image;
xlabel('Length (m)')
ylabel('Width (m)')
simple_figure()
% ------------------------------------------------------------------------------
figure;
subplot(1,2,1);
hold on;
plot(x(1:3),x(4:6),'k.','markersize',30)
plot(0,0,'k.','markersize',30)
hold off;
grid on;
axis image;
xlabel('Length (m)')
ylabel('Width (m)')
simple_figure()

subplot(1,2,2);
loglog(objfnc,'r.-','markersize',20)
% plot(objfnc,'r.-','markersize',20)
axis tight;
xlabel('Iteration #')
ylabel('Objective function')
simple_figure()
% ------------------------------------------------------------------------------
