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
np = 18;
max_neigh = 7;
neigh = zeros(np,max_neigh,'uint32');

neigh(1,1:3) = [2 , 4, 7];
neigh(2,1:5) = [1 , 3, 4, 5,  8];
neigh(3,1:4) = [2 , 5, 6, 9];
neigh(4,1:5) = [1 , 2, 7, 8, 10];
neigh(5,1:4) = [2 , 3, 8, 9];
neigh(6,1:3) = [3, 9, 12];
neigh(7,1:4) = [1, 4, 10, 13];
neigh(8,1:7) = [2, 4, 5, 9, 10, 11, 14];
neigh(9,1:7) = [3, 5, 6, 8, 11, 12, 15];
neigh(10,1:6)= [4, 7, 8, 13, 14, 16];
neigh(11,1:4)= [8, 9, 14, 15];
neigh(12,1:4)= [6, 9, 15, 18];
neigh(13,1:3)= [7, 10, 16];
neigh(14,1:6)= [8, 10, 11, 15, 16, 17];
neigh(15,1:6)= [9, 11, 12, 14, 17, 18];
neigh(16,1:4)= [10, 13, 14, 17];
neigh(17,1:4)= [14, 15, 16, 18];
neigh(18,1:3)= [12, 15, 17];

xy = randn(np*2,1);
xy = 2*[0; 6; 12; 3; 9; 15; 0; 6; 12; 3; 9; 15; 0; 6; 12; 3; 9; 15; 0; 0; 0; 2.5; 2.5; 2.5; 5; 5; 5; 7.5; 7.5; 7.5; 10; 10; 10; 12.5; 12.5; 12.5];
distao = [5.9; 4.4; 4.7; 6.3; 3.1; 4.4; 4.5; 3.3; 3.4; 4.4; 4.1; 3.9; 4.5; 4.9; 3.4; 3.8; 3.7; 3.1; 3.4; 6.3; 3.4; 3.7; 3.1; 3.2; 3.6; 3.1; 3.3; 3.5; 2.7; 3.7; 3.2; 3.8; 3.5; 3.1; 6; 3.3; 3.7; 3; 4; 6.6; 6.2];
% ------------------------------------------------------------------------------
%                    🌹🌹🌹🌹🌹🌹🌹
% ------------------------------------------------------------------------------
np = 6;
max_neigh = 4;
neigh = zeros(np,max_neigh,'uint32');

neigh(1,1:4) = [2 3 5 6];
neigh(2,1:2) = [1, 3];
neigh(3,1:4) = [1, 2, 4, 5];
neigh(4,1:2) = [3, 5];
neigh(5,1:4) = [1, 3, 4, 6];
neigh(6,1:2) = [1, 5];

xy = zeros(np*2,1);
radii=6;
for ixy=1:np
xy(ixy) = radii*cos(ixy*(2*pi/np));
xy(np+ixy) = radii*sin(ixy*(2*pi/np));
end
distao = [7.66; 8.57; 3.69; 2.99; 1.7; 7.1; 7.1; 4.815; 1.26];
% ------------------------------------------------------------------------------
[edges,nedges] = neigh2edges(np,max_neigh,neigh);
% edges
% %{
% ------------------------------------------------------------------------------
xy_init = xy;
% ------------------------------------------------------------------------------
Jac = jac_dista_(np,nedges,max_neigh,neigh,xy);
% ------------------------------------------------------------------------------
% fwd
dista = fwd_dista_(np,nedges,xy,edges);
% ⚠️ very important to do this:
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
  axis ij;
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

  % gra_(1)=0;
  % gra_(np+1)=0;

  % step size 🚶
  step_= step_dista_(np,nedges,xy,edges,distao,gra_,iter,OBJ);
  gra_ = step_*gra_;

  % update
  xy   = xy - gra_;
  iter = iter+1;

  % re-center
  xy(1:np) = xy(1:np) - xy(1);
  xy((np+1):2*np) = xy((np+1):2*np) - xy(np+1);
end
% ------------------------------------------------------------------------------
%                                  🎨
% ------------------------------------------------------------------------------
figure;
hold on;
% plot(xy_init(1:np),xy_init((np+1:2*np)),'b.','markersize',35)
plot(xy(1:np),xy((np+1:2*np)),'k.','markersize',30)
hold off;
grid on;
axis image;
axis ij;
xlabel({'Length (m)';'Parking'})
ylabel({'Warehouse';'Width (m)'})
simple_figure()

figure;
% loglog(objfnc,'r.-','markersize',20)
semilogy(10.^objfnc,'r.-','markersize',20)
axis tight;
xlabel('Iteration #')
ylabel('Objective function')
simple_figure()
% ------------------------------------------------------------------------------
%}
