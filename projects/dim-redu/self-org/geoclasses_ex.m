clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
addpath('../../opti/dense-circ/src/')
% ------------------------------------------------------------------------------
path_data = '../../../data/self-org/';
name_data = 'geoclasses.txt';
% ------------------------------------------------------------------------------
%     classess
% .---------------.
% |               |
% |               |
% |               | skills
% |               |
% .---------------.
%
% columns:
% 1. sed/strat
% 2. hydrology
% 3. geomorphology
% 4. structure
% 5. geochemistry
% 6. field
% 7. petrology
% 8. mineralogy
% 9. historical
% 10. paleontology
% 11. environment
% 12. geophysics
% 13. climate
% 14. earth materials
% 15. geology
% ------------------------------------------------------------------------------
% rows
% 1. Geologic reasoning & synthesis (dist. Obs/ interp)
% 2. Work as part of a team
% 3. Quantitative skills (algebra)
% 4. Temporal thinking
% 5. Apply skills in new scenarios (use skills from previous)
% 6. Evaluation of literature (read primary literature)
% 7. Spatial thinking (3D spatial thinking)
% 8. Written communication 
% 9. Geologic reasoning & synthesis (desc. Evidence in support argument)
% 10. Manage uncertainty (address uncertainty? when interp data)
% 11. Field Skills (make field observations)
% 12. Apply skills in new scenarios (integrate info from diff sources)
% 13. Apply skills in new scenarios (bring together geo & other knowledge)
% 14. Systems thinking (describe systems parts relationships)
% 15. Data collection & interpretation (collect & analyze data)
% 16. Spatial thinking (work with geospatial data)
% 17. Systems thinking (coplexity of scale & interactions)
% 18. Evaluation of data quality (evaluate assumptions)
% 19. Quantitative skills (statistics)
% 20. Understand societal relevance (make connections from course to lives)
% 21. Systems thinking (change has multiple effects)
% 22. Oral communication
% 23. Understand societal relevance (problem of national interest)
% 24. Evaluation of data quality (distinctions among data sources)
% 25. Systems thinking (implications & predictions)
% 26. Quantitative skills (calculus)
% 27. Field Skills (make a geologic map)
% 28. Systems thinking (current process vs. prior history)
% 29. Systems thinking (analyze feedback loops)
% 30. Systems thinking (build predictive models)
% 31. Understand societal relevance (problem of local interest)
% 32. Systems thinking (make causal maps)
% 33. Systems thinking (models of systems behavior)
% 34. Understand societal relevance (environmental justice)
% 35. Understand societal relevance (community-inspired project)
% ------------------------------------------------------------------------------
% load data
load(strcat(path_data,name_data));
d = geoclasses;
clear geoclasses;
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
[n_atributes,nd] = size(d);
d = d-repmat(mean(d),n_atributes,1);
d = normc(d);
% ------------------------------------------------------------------------------
fprintf('\nbuilding graph\n')
% build graph.
% a grid for now,
nx =  15*nd; % n_atributes*2; % 4*nd;
ny =  nx;
g = graph_grid(ny,nx);
% g = graph_torus(ny,nx);
% ------------------------------------------------------------------------------
% initialize weights for g
% --- pca
[U,S,V] = svd(d);
u = U(:,1); 
v = U(:,1) + U(:,2);
u=normc(u); 
v=normc(v);
w = w_in_grid(v,u,ny,nx);
% ------------------------------------------------------------------------------
% build self organizing map
n_iter = 30*nd;
fprintf('\nbuilding self organizing map')
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
nklust = 8;
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
