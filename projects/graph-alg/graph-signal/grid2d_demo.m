clear
close all
clc
% ------------------------------------------------------------------------------
addpath('../mesher/src/');
addpath('src/');
% ------------------------------------------------------------------------------
%                                 😷
% ------------------------------------------------------------------------------
% -- setup a simple example
% % crooked with a hole
% mask2d_ = [0 0 0 1 1 0 0; 0 1 1 1 1 0 0; 1 1 1 1 1 1 1; 1 1 0 0 0 1 1];
% crooked but hole is bigger
mask2d_ = [0 0 0 1 1 0 0; 0 1 1 1 1 0 0; 1 1 1 1 1 1 1; 1 1 0 0 0 1 1; 1 1 1 1 1 1 1];
% % full with no holes
% mask2d_ = [1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];
% mask2d_=ones(50,50);
% ------------------------------------------------------------------------------
[nz,nx]=size(mask2d_);
% ------------------------------------------------------------------------------
mask2d = zeros(nx*nz,1,'uint32');

% one for loop for three levels
for izx = 1:nx*nz
  % get x,z coordinate
  [ix,iz] = get_izx(izx,nx,nz);
  mask2d(izx)= mask2d_(iz,ix);
end
% ------------------------------------------------------------------------------
%                              🔲 & 🍇
% ------------------------------------------------------------------------------
% ◼ build structures for handling 🔲 & 🍇
tic;
n_g2m = n_g2m_(mask2d,nx,nz);
[graph2mesh,mesh2graph] = g2m_m2g(mask2d,nx,nz,n_g2m);
neigh_mesh  = neigh_mesh_(mask2d,nx,nz,n_g2m,graph2mesh);
neigh_graph = neigh_graph_(neigh_mesh,mesh2graph,n_g2m);
neigh_type  = neigh_type_(mask2d,nx,nz,n_g2m,graph2mesh);
[n_ij,n_IJ] = nIJ(n_g2m,neigh_type);
[I,J]       = IJ_(n_g2m,n_ij,n_IJ,neigh_graph);
toc;
fprintf(' just finished building all 🔲 & 🍇 structs\n\n')
% ------------------------------------------------------------------------------
%                             📶 in 🍇
% ------------------------------------------------------------------------------
x = (1:nx).';
z = (1:nz).';

circ_c = [3.5 , 2.5];
circ_r = 1;
circ_v = 2;
% ------------------------------------------------------------------------------
signal = ones(n_g2m,1);
% one for loop for two levels
for i_g2m=1:n_g2m
  % get x,y,z coordinate
  izx = graph2mesh(i_g2m);
  [ix,iz] = get_izx(izx,nx,nz);
  % ----------------------------------------------------------------------------
  % here you fill the 2d matrix.
  %
  % you can make it look like anything,
  % because you can access x, and z directly.
  %
  % below is an example on how to do that.
  % ----------------------------------------------------------------------------
  val = 1;
  % sphere filling
  radius_ = sqrt((x(ix) - circ_c(1))^2 + (z(iz) - circ_c(2))^2);
  if radius_ <= circ_r
    val = circ_v;
  end
  % ----------------------------------------------------------------------------
  % write values to 3d column
  signal(i_g2m)  = val;
end
% ------------------------------------------------------------------------------
%                                 🔵
% ------------------------------------------------------------------------------
W = ones(n_IJ,1);
V = graph_laplacian(n_g2m,n_ij,n_IJ,I,J,W);
L = sparse(I,J,V);

tic;
[E,D] = eig(full(L));
toc;
fprintf('  just finished the eigenvector decomposition\n\n')
% ------------------------------------------------------------------------------
% 🚥
filter_=zeros(n_g2m,1);

% % this filter is *exactly* a normal low-pass fourier filter,
% % but i cant figure out how it is defined 😞. see,
% % https://balcilar.medium.com/struggling-signals-from-graph-34674e699df8
% for i_g2m=1:n_g2m
%   ix=E(i_g2m,2);
%   iz=E(i_g2m,3);
%   filter_(i_g2m,1)=flt( (n_g2m+1)/2+iz , (n_g2m+1)/2+ix );
% end

d = diag(D);
filter_ = exp(-2*d/d(n_g2m));
% ------------------------------------------------------------------------------
% 🍇 𝓕
signal_= E'*signal;
% 🚥
signal_filt=filter_.*signal_;
% 🍇 𝓕^-1
signal_filt= E*signal_filt;
% ------------------------------------------------------------------------------
%                                 🎨
% ------------------------------------------------------------------------------
fprintf(' im plotting now 😴...\n\n')
% ------------------------------------------------------------------------------
% this is only for easy reference:
mask2d_ind = 1:(nx*nz);
mask2d_ind = reshape(mask2d_ind,[nz,nx]);
% this is just for visualizing
graph_ind = nan(nz,nx);
graph_ind(graph2mesh) = (1:n_g2m);
% pcolor does not plot the edges for some weird reason I do not understand
graph_ind = [graph_ind ; nan(1,nx)];
graph_ind = [graph_ind , nan(nz+1,1)];

% plot the 📶 in 🍇
signal_plot  = nan(nz,nx);
signal_plot_ = nan(nz,nx);
signal_plot_f= nan(nz,nx);
ieig = 2;
E_vec = nan(nz,nx);
filter_plot = nan(nz,nx);
% one for loop for two levels
for i_g2m=1:n_g2m
  % get x,y,z coordinate
  izx = graph2mesh(i_g2m);
  [ix,iz] = get_izx(izx,nx,nz);
  signal_plot(iz,ix)  = signal(i_g2m);
  signal_plot_(iz,ix) = signal_(i_g2m);
  signal_plot_f(iz,ix)= signal_filt(i_g2m);
  E_vec(iz,ix)= E(i_g2m,ieig);
  filter_plot(iz,ix)= filter_(i_g2m);
end
% pcolor does not plot the edges for some weird reason I do not understand
signal_plot  = [signal_plot ; nan(1,nx)];
signal_plot  = [signal_plot , nan(nz+1,1)];
signal_plot_ = [signal_plot_ ; nan(1,nx)];
signal_plot_ = [signal_plot_ , nan(nz+1,1)];
signal_plot_f= [signal_plot_f ; nan(1,nx)];
signal_plot_f= [signal_plot_f , nan(nz+1,1)];
E_vec= [E_vec ; nan(1,nx)];
E_vec= [E_vec , nan(nz+1,1)];
filter_plot= [filter_plot ; nan(1,nx)];
filter_plot= [filter_plot , nan(nz+1,1)];
% ------------------------------------------------------------------------------
figure;
subplot(131)
fancy_imagesc(mask2d_)
title('Grid mask')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

subplot(132)
fancy_imagesc(mask2d_ind)
colormap(rainbow2_cb(1))
title('Grid indexes')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

subplot(133)
fancy_pcolor(graph_ind)
colormap(rainbow2_cb(1))
title('Graph indexes')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
figure;
subplot(121)
fancy_imagesc(full(L))
colormap(rainbow2_cb(1))
title('Laplacian matrix')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

subplot(122)
fancy_pcolor(E_vec)
caxis( [min(E(:)) , max(E(:)) ] )
colormap(rainbow2_cb(1))
title(strcat('eigenvector #',num2str(ieig)))
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure
% ------------------------------------------------------------------------------
figure;
subplot(221)
fancy_pcolor(signal_plot)
colormap(rainbow2_cb(1))
caxis([min(signal) max(signal)])
title('Signal')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure

subplot(222)
fancy_pcolor(signal_plot_f)
colormap(rainbow2_cb(1))
caxis([min(signal) max(signal)])
title('Filtered signal')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

subplot(223)
fancy_pcolor(signal_plot_)
colormap(rainbow2_cb(1))
title('Spectral signal')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

subplot(224)
fancy_pcolor(filter_plot)
colormap(rainbow2_cb(1))
title('Filter')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
