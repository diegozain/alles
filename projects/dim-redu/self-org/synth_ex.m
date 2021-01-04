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
subplot(121)
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
subplot(122)
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
