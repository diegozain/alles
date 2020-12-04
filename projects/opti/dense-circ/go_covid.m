% ------------------------------------------------------------------------------
% 
% compute 2D distancing due to covid using a modified nbody solver
%
%  http://hydra.nat.uni-magdeburg.de/packing/cci/
% ------------------------------------------------------------------------------
clear
close all
clc
% ..............................................................................
addpath('src/')
% ------------------------------------------------------------------------------
% choosing the distance for separation is not straight forward because
% at large iterations the gradient is weak and wont update no more.
%
% nr=5  -> d=3  -> min distance=2
% nr=15 -> d=6  -> min distance=2
% nr=20 -> d=7  -> min distance=2
% nr=50 -> d=13 -> min distance=2
% nr=100-> d=22 -> min distance=2
% 
% a linear fit for min distance 2 is:
% 
% d = 0.195*nr + 2.785;
% 
% although at small nr it resembles more a square root,
% so this gives a lower estimate!!
nr= 19;
d = 0.195*nr + 2.785 + 0.5; % (m)
% ------------------------------------------------------------------------------
% total number of iterations
niter= 500; 
% for plotting the path of the solution: every niter__ plot a solution
niter__=2;
% ------------------------------------------------------------------------------
fprintf('\n\n   %i nodes and min distance %2.2d  \n\n',nr,d);
% ------------------------------------------------------------------------------
rng(3);
r = zeros(nr,2);
r = randn(nr,2) * (d/5);
r_= r;
% ------------------------------------------------------------------------------
% gradient descent
[r,R,E,steps_] = covid_neigh(r,d,niter);
% ------------------------------------------------------------------------------
% % hessian 
% r = covid_neigh_(r,d);
% ------------------------------------------------------------------------------
m = geome_median(r);
D = dist_mat(r);
fprintf('\n max and min distances: %2.2d , %2.2d\n\n',max(D(:)),min( D(D(:)>0) ))
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(D)
colormap(rainbow())
ylabel('Point #')
xlabel('Point #')
title('Distance Matrix (m)')
simple_figure()

figure;
hold on
plot(r(:,1),r(:,2),'r.','markersize',50)
plot(r_(:,1),r_(:,2),'k.','markersize',20)
plot(m(1),m(2),'b.','markersize',30)
hold off
axis tight
axis image
xlabel('Length (m)')
ylabel('Width (m)')
title('Original vs Optimal')
simple_figure()

figure;
plot(1:size(E,1),E/nr,'.-','markersize',10)
axis tight
xlabel('Iteration #')
ylabel('Objective functions per node')
simple_figure()

figure;
plot(1:size(E,1),log10(sum(E,2)),'.-','markersize',10)
axis tight
xlabel('Iteration #')
ylabel('Objective function')
simple_figure()

figure;
plot(1:size(steps_,1),steps_,'.-','markersize',10)
axis tight
xlabel('Iteration #')
ylabel('Step sizes per node')
simple_figure()
% ------------------------------------------------------------------------------
%%{
niter = size(R,3);
niter_= floor(niter/niter__);

xmax=max(R(:,1,:));
xmax=max(xmax);
xmin=min(R(:,1,:));
xmin=min(xmin);
ymax=max(R(:,2,:));
ymax=max(ymax);
ymin=min(R(:,2,:));
ymin=min(ymin);

figure;
hold on
plot(R(:,1,1),R(:,2,1),'k.','markersize',20)
for i_=1:niter_
 plot(R(:,1,(i_-1)*niter__+1),R(:,2,(i_-1)*niter__+1),'.','markersize',10)
 xlim([xmin xmax])
 ylim([ymin ymax])
 xlabel('Length (m)')
 ylabel('Width (m)')
 title('Original vs Optimal')
 simple_figure()
 % pause;
end
plot(R(:,1,niter),R(:,2,niter),'r.','markersize',40)
hold off;
%}
% ------------------------------------------------------------------------------
r(:,1) = r(:,1)-repmat(m(1),nr,1);
r(:,2) = r(:,2)-repmat(m(2),nr,1);

figure;
hold on
plot(r(:,1),r(:,2),'r.','markersize',50)
plot(0,0,'b.','markersize',30)
hold off
axis tight
axis image
xlabel('Length (m)')
ylabel('Width (m)')
title('Optimal')
simple_figure()
% ------------------------------------------------------------------------------
