clear
close all
clc
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
np = 4;
max_neigh = 3;
% ------------------------------------------------------------------------------
neigh = zeros(np,max_neigh,'uint32');

neigh(1,1:3) = [2 , 3, 4];
neigh(2,1:3) = [1 , 3, 4];
neigh(3,1:3) = [1 , 2, 4];
neigh(4,1:3) = [1 , 2, 3];
% ------------------------------------------------------------------------------
[edges,nedges] = neigh2edges(np,max_neigh,neigh);
% edges
% %{
% ------------------------------------------------------------------------------
xy = zeros(np*2,1);
xy = [4.8; 3; -1; 0; 0; 4; 4; 0];
% ------------------------------------------------------------------------------
xy_init = xy;
% ------------------------------------------------------------------------------
Jac = jac_dista_(np,nedges,max_neigh,neigh,xy);
% ------------------------------------------------------------------------------
% fwd
dista = fwd_dista_(np,nedges,xy,edges);

distao = [4.5; 6.5;  4.8; 3.0; 3.6; 3.85];
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

  gra_(np)=0;
  gra_(np*2)=0;

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
xlabel('Iteration #')
ylabel('Objective function')
simple_figure()
% ------------------------------------------------------------------------------
%}
