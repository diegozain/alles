clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% total number of electrodes
nelectrodes  = 48; % 64 48
% ------------------------------------------------------------------------------
fprintf('total # of electrodes = %i\n',nelectrodes)
% ------------------------------------------------------------------------------
ielect_start = 4;
ielect_inter = 4;
% electrode indexes that will be Tx
Tx_ = uint32((ielect_start:ielect_inter:(nelectrodes*0.5)).');
% electrode indexes that will be Rx
Rx_ = uint32(((nelectrodes*0.5 + 1):nelectrodes).');
% ------------------------------------------------------------------------------
% get all abmn pairs
abmn = xbore_getall_(Tx_,Rx_);
nabmn= size(abmn,1);
% ------------------------------------------------------------------------------
% declare all electrode positions
xmin = 1;
xmax = 4;
zmin = 1;
zmax = 1 * nelectrodes*0.5;
electrodes=[xmin*ones(nelectrodes/2,1) (linspace(zmin,zmax,nelectrodes/2)).'; xmax*ones(nelectrodes/2,1) (linspace(zmin,zmax,nelectrodes/2)).'];
% ------------------------------------------------------------------------------
fprintf('total # of abmn quadruples = %i\n',nabmn)
% ------------------------------------------------------------------------------
% get all pseudo locations
tic;
pseud= xbore_pseudo(electrodes,abmn);
toc;
% ------------------------------------------------------------------------------
fprintf(' -- just finished all pseudo locations -- \n\n')
% ------------------------------------------------------------------------------
figure;
subplot(1,3,1)
hold on;

plot(electrodes(Rx_(:),1),electrodes(Rx_(:),2),'k.','markersize',40)
plot(electrodes(Tx_(:),1),electrodes(Tx_(:),2),'.','markersize',40,'color',[0.5,0.5,0.5])

plot(pseud(:,1),pseud(:,2),'.','color',[0.9098,0.2196,0.1451],'markersize',20)

hold off;
axis ij;
axis image;
xlim([xmin-2,xmax+2])
ylim([0,zmax+2])
xlabel('Length (m)')
ylabel('Depth (m)')
title('All AB.MN pseudo locations')
simple_figure()
% ------------------------------------------------------------------------------
% 1. get repeated elements in pseud,
% 2. count # of clusters,
% 3. keep track of what clusters are,
% 4. build clusters.
% 
% #3 proceeds to make a vector 'klusters' where entries are how many elements 
% repeated elements belong to a cluster (they are ordered that way).
% 
% example (not real example):
% 
% repeated    = [2; 4; 6; 10; 11; 15]
% if clusters = [2 4] , [6 10 11 15], then
%    klusters = [2; 4]
% ------------------------------------------------------------------------------
tic;
[klusters_,repeated] = xbore_clusters(abmn,pseud);
toc;
% ------------------------------------------------------------------------------
nrepeat = numel(repeated);
nklu = size(klusters_,1);
% ------------------------------------------------------------------------------
fprintf('total # of repeated abmn   = %i\n',nrepeat)
fprintf('repeated / total           = %2.2f\n\n',nrepeat / nabmn)
fprintf('total # of clusters        = %i\n\n',nklu)
% ------------------------------------------------------------------------------
subplot(1,3,2)
hold on;
plot(electrodes(Rx_(:),1),electrodes(Rx_(:),2),'k.','markersize',40)
plot(electrodes(Tx_(:),1),electrodes(Tx_(:),2),'.','markersize',40,'color',[0.5,0.5,0.5])
% iabmn=1;
for iklu=1:nklu 
  % fprintf('cluster #%i\n',iklu)
  nklu_=size(klusters_{iklu},1);
  
  % this part is to plot each cluster with a different color and a different
  % size.
  % I dont really know what determines the size of a cluster.
  % the larger the # of electrodes, the more cluster configurations you get,
  % and not just in the obvious case of repeating the previous configurations, 
  % you also get new configurations! 
  % and these new configs are sometimes of smaller size than the previously
  % found.
  % 
  % so far, the clusters are of size:
  % 3, 5, 6, 9, 10, 14, 15, 19, 20, 21, 28, 36, 45, 50, 55, 66, 78
  msize = nklu_*3;
  colo=[0.0941,0.3961,0.8667];

  plot(pseud(klusters_{iklu},1),pseud(klusters_{iklu},2),'.','markersize',msize,'color',colo)

  % for iklu_=1:nklu_
  %   elec=electrodes(abmn(klusters_{iklu}(iklu_),:),:);
  %   plot(elec(:,1),elec(:,2),'k.','markersize',40)
  %   % pause;
  %   plot(elec(:,1),elec(:,2),'.','markersize',40,'color',[0.5,0.5,0.5])
  % end
end
hold off;
axis ij;
axis image;
xlim([xmin-2,xmax+2])
ylim([0,zmax+2])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Repeated locations')
simple_figure()
% ------------------------------------------------------------------------------
% 
% 
% visualize the data
% 
% 
% ------------------------------------------------------------------------------
subplot(1,3,3)
hold on;
% plot electrodes
plot(electrodes(Rx_(:),1),electrodes(Rx_(:),2),'k.','markersize',40)
plot(electrodes(Tx_(:),1),electrodes(Tx_(:),2),'.','markersize',40,'color',[0.5,0.5,0.5])
% plot all pseudo positions (including repetitions)
scatter(pseud(:,1),pseud(:,2),70,'r','filled')
% plot only repeated positions (clusters)
for iklu=1:nklu 
  % fprintf('cluster #%i\n',iklu)
  nklu_=size(klusters_{iklu},1);
  % first, let's paint white all the repeated locs we already plotted
  scatter(pseud(klusters_{iklu},1),pseud(klusters_{iklu},2),70,'w','filled')
  % ----------------------------------------------------------------------------
  path_circs   = strcat('../../../data/cci_coords/cci',num2str(nklu_),'.txt');
  dense_circle = load(path_circs);
  % the guy who wrote these put an index up front, so we have to remove that
  dense_circle = dense_circle(:,2:3);
  
  pseud_ = pseud(klusters_{iklu},:) + 0.1*dense_circle;
  
  scatter(pseud_(:,1),pseud_(:,2),20,'b','filled')
  % ----------------------------------------------------------------------------
end
hold off;
axis ij;
axis image;
xlim([xmin-2,xmax+2])
ylim([0,zmax+2])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Pseudo sections')
simple_figure()
% ------------------------------------------------------------------------------
% 
% 
% visualize just one abmn pseudo location and its sensitivity
% 
% 
% ------------------------------------------------------------------------------
nx = 6e2;
nz = 6e2;

x=linspace(xmin-2,xmax+2,nx);
z=linspace(0,zmax+2,nz);

[X,Z] = meshgrid(x,z);

dx=x(2)-x(1);
dz=z(2)-z(1);
% ------------------------------------------------------------------------------
figure;
for iexample=1:16
  % choose abmn configuration
  iabmn=randi(nabmn);

  source_p = electrodes(abmn(iabmn,1),:);
  source_n = electrodes(abmn(iabmn,2),:);

  rec_p = electrodes(abmn(iabmn,3),:);
  rec_n = electrodes(abmn(iabmn,4),:);
  % ----------------------------------------------------------------------------
  psi_  = sensitivity_3dDC(X,Z,dx,dz,source_p,source_n,rec_p,rec_n);
  psi_  = normali(psi_);
  % ----------------------------------------------------------------------------
  subplot(2,8,iexample)
  fancy_imagesc(psi_,x,z)
  colormap(rainbow2_cb(1))
  caxis(5e-4*[-1 1])
  colorbar('off')
  hold on;
  plot(pseud(iabmn,1),pseud(iabmn,2),'.','markersize',40,'color',[0.1529,0.7686,0.0980])
  plot(pseud(iabmn,1),pseud(iabmn,2),'w.','markersize',20)
  hold off;
  set(gca,'xtick',[])
  set(gca,'ytick',[])
end
simple_figure()
% ------------------------------------------------------------------------------
