clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% data
% dimensions of (attributes x data points)
d = [255 0 0 0; 0 255 0 20; 0 0 255 200];
[n_atributes,nd] = size(d);
% ------------------------------------------------------------------------------
figure;
subplot(131)
fancy_imagesc(d)
colormap(rainbow2(1))
colorbar('off')
title('Data')
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
d = d-repmat(mean(d),n_atributes,1);
d = normc(d);
% ------------------------------------------------------------------------------
% build graph.
% a grid for now,
nx =  3*nd;
ny =  nx;
g = graph_grid(ny,nx);
% g = graph_torus(ny,nx);
% ------------------------------------------------------------------------------
% initialize weights for g
% --- pca
[U,S,V] = svd(d);
u = U(:,1); 
v = U(:,2);
u=normc(u); 
v=normc(v);
w = w_in_grid(v,u,ny,nx);
% ------------------------------------------------------------------------------
% build self organizing map
n_iter = 30*nd;
fprintf('\nbuilding self organizing')
fprintf('\n                          size: %i by %i pixels',ny,nx)
fprintf('\n                          iter: %i\n\n',n_iter)
[w,u_mat,d_in_g] = selforgmapi(d,g,w,n_iter);
% ------------------------------------------------------------------------------
% return graph and data to grid
u_mat = reshape(u_mat,ny,nx);
d_in_g_ = zeros(nd,2);
for i_=1:nd
  [a,b]=ind2sub([ny,nx],d_in_g(i_));
  d_in_g_(i_,:) = [a,b];
end
% ------------------------------------------------------------------------------
% see
subplot(132)
imagesc(u_mat);
colormap(flipud(bone)); 
caxis([0 1])
hold on
plot( d_in_g_(:,2),d_in_g_(:,1),'r.','markersize',80 )
a = [1:nd]'; b = num2str(a); cc = cellstr(b);
text(d_in_g_(:,2), d_in_g_(:,1), cc,'FontSize',14,'HorizontalAlignment','center');
hold off
axis square
ax=gca;
ax.XTick=[];
ax.YTick=[];
title('SOM of data')
simple_figure()
% ------------------------------------------------------------------------------
% 
% cool cluster plot
% 
% 
% ------------------------------------------------------------------------------
fprintf('\n cool plotting with clusters now\n')
nklust = 3;
iklust = kmeans(d_in_g_,nklust);
% ------------------------------------------------------------------------------
% get # of elements for each cluster
nik=zeros(nklust,1);
for i_=1:nklust
  nik(i_) = numel(iklust(iklust==i_));
end

% get elements ID for each cluster
clusters = cell(nklust,1);
for i_=1:nklust
  clusters{i_} = find(iklust==i_);
end
% ------------------------------------------------------------------------------
rng(3);
nr= nklust;
niter= 200;
dista=4;
% ------------------------------------------------------------------------------
% cluster centers in xz coordinates for plotting
r = zeros(nr,2);
r = randn(nr,2) * (dista/15);
[r,~,~,~] = covid_neigh(r,dista,500);
m = geome_median(r);
r(:,1) = r(:,1)-repmat(m(1),nr,1);
r(:,2) = r(:,2)-repmat(m(2),nr,1);
% ------------------------------------------------------------------------------
% for each cluster and for each point in that cluster get xz coordinates 
R_=cell(nklust,1);
for i_=1:nklust
  r_ = zeros(nik(i_),2);
  r_ = randn(nik(i_),2) * (dista/5);
  [r_,~,~,~] = covid_neigh(r_,dista/(max(nik)),niter);
  m = geome_median(r_);
  r_(:,1) = r_(:,1)-repmat(m(1),nik(i_),1) + r(i_,1);
  r_(:,2) = r_(:,2)-repmat(m(2),nik(i_),1) + r(i_,2);
  
  R_{i_} = r_;
end
% ------------------------------------------------------------------------------
subplot(133)
hold on
% % centers of clusters
% plot(r(:,1),r(:,2),'b.','markersize',50)
for i_=1:nklust
  % clusters
  plot(R_{i_}(:,1),R_{i_}(:,2),'r.-','markersize',80,'linewidth',5)
  % annotations
  a = clusters{i_};
  b = num2str(a); 
  cc = cellstr(b);
  text(R_{i_}(:,1),R_{i_}(:,2), cc,'FontSize',14,'HorizontalAlignment','center');
end
hold off
axis square
ax=gca;
ax.XTick=[];
ax.YTick=[];
title('Clusters of SOM')
simple_figure()
% ------------------------------------------------------------------------------