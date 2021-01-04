clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
path_data = '../../../data/self-org/';
name_data = 'zodiac.txt';
% ------------------------------------------------------------------------------% below is the data. 
% each row is an attribute and each column a person.
% 1.Nayani 
% 2.Hamid 
% 3.Zach 
% 4.Monica 
% 5.Thomas 
% 6.Nick 
% 7.Anna 
% 8.Jake 
% 9.Silvia 
% 10.Karun 
% 11.Josh 
% 12.Anand 
% 13.Diego 
% 14.Nancy
%
% indi 5 0 4 0 3 0 1 0 0 0 0 0 0 0
% focu 0 5 0 0 0 0 0 0 1 0 0 0 0 0
% rest 0 0 5 4 0 0 0 1 0 3 0 1 0 0
% inpu 0 0 0 5 0 0 0 0 0 1 0 4 0 0
% conn 0 0 3 0 5 0 0 0 0 0 0 0 0 0
% inte 0 0 0 0 0 5 0 0 0 0 2 2 0 0
% lear 0 0 0 2 4 3 5 0 0 5 4 5 5 0
% anal 0 0 1 0 0 0 0 5 0 0 0 0 0 1
% posi 0 0 0 0 0 0 0 0 5 0 0 0 0 0
% rela 0 0 0 1 0 2 0 0 0 0 5 0 0 0
% achi 4 0 0 0 2 0 0 0 3 0 0 0 0 3
% stra 0 4 0 0 0 0 0 0 0 0 0 0 2 0
% deli 3 0 0 0 0 4 0 0 0 0 0 0 0 0
% harm 0 0 0 0 0 0 4 0 0 4 0 0 0 0
% idea 0 2 0 0 0 0 0 4 0 0 0 0 0 0
% sign 0 1 0 0 0 0 0 0 4 0 0 0 0 4
% acti 0 0 0 0 0 0 0 0 0 0 0 0 4 0
% futu 0 3 0 0 0 0 0 0 0 0 0 0 1 2
% comm 0 0 0 3 0 0 0 0 0 0 0 0 0 0
% resp 1 0 0 0 1 1 3 0 2 0 1 3 0 5
% cons 0 0 0 0 0 0 2 3 0 0 0 0 3 0
% adap 2 0 0 0 0 0 0 0 0 0 3 0 0 0
% arra 0 0 2 0 0 0 0 0 0 0 0 0 0 0
% comp 0 0 0 0 0 0 0 2 0 0 0 0 0 0
% empa 0 0 0 0 0 0 0 0 0 2 0 0 0 0
% 
% the data is in txt format in data.txt
% ------------------------------------------------------------------------------
load(strcat(path_data,name_data));
d = zodiac; 
clear data;
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
[n_atributes,nd] = size(d);
d = d-repmat(mean(d),n_atributes,1);
% d = d / max(abs(d(:)));
d = normc(d);
% ------------------------------------------------------------------------------
% build graph.
% a grid for now,
nx =  3*nd; % n_atributes*2; % 4*nd;
ny =  nx;
g = graph_grid(ny,nx);
% g = graph_torus(ny,nx);
% ------------------------------------------------------------------------------
% initialize weights for g
% --- pca fancy
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