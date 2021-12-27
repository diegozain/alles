clear
close all
clc
% ------------------------------------------------------------------------------
addpath('../mesher/src/');
addpath('src/');
% ------------------------------------------------------------------------------
%                                🍇
% ------------------------------------------------------------------------------
ngraph = 100;
nei_max = 2; % because the graph is cyclic
% ------------------------------------------------------------------------------
neigh_graph = neigh_graph_cyclic_(ngraph);
[n_ij,n_IJ] = nIJ_graph(ngraph,nei_max,neigh_graph);
[I,J] = IJ_graph(ngraph,nei_max,n_ij,n_IJ,neigh_graph);
% ------------------------------------------------------------------------------
%                             📶 in 🍇
% ------------------------------------------------------------------------------
signal = randn([ngraph,1]);

% signal = zeros(ngraph,1);
% signal(fix(ngraph/3)) = 1;
% signal(fix(2*ngraph/3)) = 1;
% signal(fix(ngraph)) = 1;
% ------------------------------------------------------------------------------
%                                 🔵
% ------------------------------------------------------------------------------
W = ones(n_IJ,1);
V = graph_laplacian(ngraph,n_ij,n_IJ,I,J,W);
L = sparse(I,J,V);
% ------------------------------------------------------------------------------
tic;
[E,D] = eig(full(L));
toc;
fprintf('  just finished the eigenvector decomposition\n\n')
% ------------------------------------------------------------------------------
% 🚥
filter_=zeros(ngraph,1);

% % this filter is *exactly* a normal low-pass fourier filter,
% % but i cant figure out how it is defined 😞. see,
% % https://balcilar.medium.com/struggling-signals-from-graph-34674e699df8
% for i_g2m=1:ngraph
%   ix=E(i_g2m,2);
%   iz=E(i_g2m,3);
%   filter_(i_g2m,1)=flt( (ngraph+1)/2+iz , (ngraph+1)/2+ix );
% end

d = diag(D);
filter_ = exp(-2*d/d(ngraph));
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
% eigenvector of laplacian to plot
ieig = 2;
% thickness of nodes in plot
thic = 50;
% ------------------------------------------------------------------------------
V = graph_adjacent(ngraph,n_ij,n_IJ,I,J,W);
A = sparse(I,J,V);
% ------------------------------------------------------------------------------
thet = linspace(0,2*pi,ngraph+1);
thet = thet.';
thet(1) = [];
nodes = [cos(thet) , sin(thet)];
% ------------------------------------------------------------------------------
% 🍇
plot_graph(nodes,A,thic/2);
% ------------------------------------------------------------------------------
%  📶 in 🍇
plot_graph_(nodes,A,signal,thic);
% colormap(rainbow())
caxis([min(signal),max(signal)])
title('Raw')

plot_graph_(nodes,A,signal_filt,thic);
% colormap(rainbow())
caxis([min(signal),max(signal)])
title('Filtered')

plot_graph_(nodes,A,signal_,thic);
title('Spectral')

plot_graph_(nodes,A,filter_,thic);
title('Filter')

plot_graph_(nodes,A,E(:,ieig),thic);
caxis( [min(E(:)) , max(E(:)) ] )
title(strcat('eigenvector #',num2str(ieig)))

plot_graph_(nodes,A,d,thic);
title('Eigenvalues')
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(full(L))
colormap(rainbow2_cb(1))
title('Laplacian matrix')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
