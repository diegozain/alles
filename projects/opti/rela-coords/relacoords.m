clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src/')
% ------------------------------------------------------------------------------
% np        = # of points
% max_neigh = maximum # of neighbors a node may have
% neigh     = neighbor list of size (np × max_neigh) ℕ
% edges     = numbered edge collection of size (# of edges × 2) ℕ
%             a construct of neigh
% xy        = xy coordinates of points,
%             of size (np*2 × 1) ℝ, first half is x, second is y
% dist      = list of distances between all points
%             with the same order as 'edges' (# of edges × 1) ℝ>0
% ------------------------------------------------------------------------------
np = 7;
max_neigh = 5;
% ------------------------------------------------------------------------------
neigh = zeros(np,max_neigh,'uint32');

neigh(1,1:5) = [2, 3, 4, 5, 6];
neigh(2,1:2) = [1, 3];
neigh(3,1:5) = [1, 2, 4, 6, 7];
neigh(4,1:4) = [1, 3, 5, 6];
neigh(5,1:3) = [1, 4, 6];
neigh(6,1:5) = [1, 3, 4, 5, 7];
neigh(7,1:4) = [3, 4, 5, 6];
% ------------------------------------------------------------------------------
[edges,nedges] = neigh2edges(np,max_neigh,neigh);
% edges
% %{
% ------------------------------------------------------------------------------
xy = zeros(np*2,1);
xy = [0; 0; 0; 15; 12.7; 18.5; 12.7; 0; 0; 0; 0; -2.7; -1.2; -4];
% ------------------------------------------------------------------------------
load('..\..\..\..\dcip\field\protocols\data\geome\borexyextend.mat');
xy(1) = borexyextend(15,1);
xy(2) = borexyextend(17,1);
xy(3) = borexyextend(18,1);
xy(8) = borexyextend(15,2);
xy(9) = borexyextend(17,2);
xy(10)= borexyextend(18,2);
% ------------------------------------------------------------------------------
xy_init = xy;
% ------------------------------------------------------------------------------
Jac = jac_dista_(np,nedges,max_neigh,neigh,xy);
% ------------------------------------------------------------------------------
% fwd
dista = fwd_dista_(np,nedges,xy,edges);
% ------------------------------------------------------------------------------
distao = [ 0; 0; 3.95; 6.15; 6.6; 0; 4.1; 4.2; 8.3; 2.35; 3.5; 3.6; 4.65];
% ------------------------------------------------------------------------------
distao(1) = sqrt((borexyextend(15,1) - borexyextend(17,1)).^2 + (borexyextend(15,2) - borexyextend(17,2)).^2);
distao(2) = sqrt((borexyextend(15,1) - borexyextend(18,1)).^2 + (borexyextend(15,2) - borexyextend(18,2)).^2);
distao(6) = sqrt((borexyextend(17,1) - borexyextend(18,1)).^2 + (borexyextend(17,2) - borexyextend(18,2)).^2);
% ------------------------------------------------------------------------------
distao = distao.^2;
% ------------------------------------------------------------------------------
%
%                              👈 inv 👈
%
% ------------------------------------------------------------------------------
niter  = 100;
OBJ = 'rms_log';
objfnc_= Inf;
objfnc = zeros(niter,1);
tol_error = -11;%1e-5;
% ------------------------------------------------------------------------------
iter   = 1;
while (objfnc_>tol_error & iter<=niter)
  % 🎨
  figure(1001);
  hold on;
  plot(xy(1:np),xy((np+1:2*np)),'.','markersize',20)
  hold off;
  axis image;
  xlabel('Length (m)')
  ylabel('Width (m)')
  simple_figure()

  % forward model 👉
  dista = fwd_dista_(np,nedges,xy,edges);

  % obj function 👇
  [objfnc_,err_]= obj_dista(dista, distao, OBJ);
  objfnc(iter)  = objfnc_;

  % bacwkard model 👈
  Jac  = jac_dista_(np,nedges,max_neigh,neigh,xy);
  gra_ = Jac * err_;

  % ⛔ these points have no update
  gra_(1)=0;
  gra_(2)=0;
  gra_(3)=0;
  gra_(np+1)=0;
  gra_(np+2)=0;
  gra_(np+3)=0;

  % step size 🚶
  step_= step_dista_(np,nedges,xy,edges,distao,gra_,iter,OBJ);
  gra_ = step_*gra_;

  % update
  xy   = xy - gra_;
  iter = iter+1;

  % % re-center
  % xy(1:np) = xy(1:np) - xy(1);
  % xy((np+1):2*np) = xy((np+1):2*np) - xy(np+1);
end
% ------------------------------------------------------------------------------
%                                  🎨
% ------------------------------------------------------------------------------
figure;
hold on;
plot(xy_init(1:np),xy_init((np+1:2*np)),'b.','markersize',35)
plot(xy(1:np),xy((np+1:2*np)),'k.','markersize',30)
hold off;
grid on;
axis image;
xlabel('Length (m)')
ylabel('Width (m)')
simple_figure()

figure;
loglog(objfnc,'r.-','markersize',20)
% plot(objfnc,'r.-','markersize',20)
axis tight;
axis square;
grid on;
xlabel('Iteration #')
ylabel('Objective function')
simple_figure()
% ------------------------------------------------------------------------------


% ------------------------------------------------------------------------------
borexyextend_=borexyextend;
borexyextend_(19:22,1) = xy(4:7);
borexyextend_(19:22,2) = xy(11:14);

figure;
plot(borexyextend_(:,1),borexyextend_(:,2),'.','markersize',35)
grid on;
axis image;
xlabel('Length (m)')
ylabel('Width (m)')
simple_figure()
% ------------------------------------------------------------------------------
%}
