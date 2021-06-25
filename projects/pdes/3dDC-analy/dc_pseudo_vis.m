clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% total number of electrodes
nelectrodes=18; % 32; 
% declare all electrode positions
electrodes=[1*ones(nelectrodes/2,1) (linspace(1,nelectrodes*0.5,nelectrodes/2)).'; 4*ones(nelectrodes/2,1) (linspace(1,nelectrodes*0.5,nelectrodes/2)).'];
% # of electrodes that will be Tx
nTx=nelectrodes/2;
% ------------------------------------------------------------------------------
% get all abmn pairs
abmn = xbore_getall(nelectrodes,nTx);
nabmn= size(abmn,1);
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
for ie=1:nelectrodes
plot(electrodes(ie,1),electrodes(ie,2),'k.','markersize',40)
end
for iabmn=1:nabmn
plot(pseud(iabmn,1),pseud(iabmn,2),'r.','markersize',20)
% plot(pseud(iabmn,1),pseud(iabmn,2),'.')
end
hold off;
axis image;
axis ij
xlim([0.75,4.25])
ylim([0,nelectrodes*0.5+1])
xlabel('Length (m)')
ylabel('Depth (m)')
title('All AB.MN pseudo locations')
simple_figure()
% ------------------------------------------------------------------------------
% this is just to make sure the algo supports abmn pairs that are not ordered.
ipermu= randperm(nabmn);
abmn  = abmn(ipermu,:);
pseud = pseud(ipermu,:);
% ------------------------------------------------------------------------------
% 1. get repeated elements in pseud,
% 2. count # of clusters,
% 3. keep track of what clusters are,
% 4. build clusters.
% 
% #3 proceeds to make a vector 'klusters' where entries are how many elements 
% repeated elements belong to a cluster (they are ordered that way).
% example (not real example):
% 
% repeated = [2; 4; 6; 10; 11; 15]
% if clusters = [2 4] , [6 10 11 15], then
% klusters = [2; 4]
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
for ie=1:nelectrodes
plot(electrodes(ie,1),electrodes(ie,2),'.','markersize',40,'color',[0.5,0.5,0.5])
end
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
  if (nklu_==3)
    colo=[0.5,0,0];
    msize=10;
  elseif (nklu_==5)
    colo=[1,0,0];
    msize=15;
  elseif (nklu_==6)
    colo=[0,0.5,0];
    msize=20;
  elseif (nklu_==9)
    colo=[0,1,0];
    msize=22;
  elseif (nklu_==10)
    colo=[0,0,0.5];
    msize=24;
  elseif (nklu_==14)
    colo=[0,0,1];
    msize=26;
  elseif (nklu_==15)
    colo=[0.5,0.5,0];
    msize=28;
  elseif (nklu_==19)
    colo=[1,0.5,0];
    msize=30;
  elseif (nklu_==20)
    colo=[1,1,0];
    msize=32;
  elseif (nklu_==21)
    colo=[1,1,0.5];
    msize=34;
  elseif (nklu_==28)
    colo=[0.5,0,0.5];
    msize=36;
  elseif (nklu_==36)
    colo=[1,0,0.5];
    msize=38;
  elseif (nklu_==45)
    colo=[1,0,1];
    msize=40;
  elseif (nklu_==50)
    colo=[1,0.5,1];
    msize=42;
  elseif (nklu_==55)
    colo=[0,0.5,0.5];
    msize=44;
  elseif (nklu_==66)
    colo=[0,1,0.5];
    msize=46;
  elseif (nklu_==78)
    colo=[0,1,1];
    msize=48;
  else
    colo=[0,0,0];
    msize=50;
    fprintf('new cluster size! =%i\n',nklu_)
  end

  plot(pseud(klusters_{iklu},1),pseud(klusters_{iklu},2),'.','markersize',msize,'color',colo)

  for iklu_=1:nklu_
    elec=electrodes(abmn(klusters_{iklu}(iklu_),:),:);
    plot(elec(:,1),elec(:,2),'k.','markersize',40)
    % pause;
    plot(elec(:,1),elec(:,2),'.','markersize',40,'color',[0.5,0.5,0.5])
  end
end
hold off;
axis ij
axis image;
xlim([0.75,4.25])
ylim([0,nelectrodes*0.5+1])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Repeated locations')
simple_figure()
% ------------------------------------------------------------------------------
subplot(1,3,3)
hold on;
% plot electrodes
scatter(electrodes(:,1),electrodes(:,2),150,'k','filled')
% plot all pseudo positions (including repetitions)
scatter(pseud(:,1),pseud(:,2),80,'r','filled')
% plot only repeated positions (clusters)
for iklu=1:nklu 
  % fprintf('cluster #%i\n',iklu)
  nklu_=size(klusters_{iklu},1);
  % first, let's paint white all the repeated locs we already plotted
  scatter(pseud(klusters_{iklu},1),pseud(klusters_{iklu},2),100,'w','filled')
  % ----------------------------------------------------------------------------
  path_circs   = strcat('../../../data/cci_coords/cci',num2str(nklu_),'.txt');
  dense_circle = load(path_circs);
  % the guy who wrote these put an index up front, so we have to remove that
  dense_circle = dense_circle(:,2:3);
  
  pseud_ = pseud(klusters_{iklu},:) + 0.05*sqrt(nklu_-2)*dense_circle;
  
  scatter(pseud_(:,1),pseud_(:,2),20,'b','filled')
  % ----------------------------------------------------------------------------
end
hold off;
axis image;
axis ij
xlim([0.75,4.25])
ylim([0,nelectrodes*0.5+1])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Pseudo sections')
simple_figure()
% ------------------------------------------------------------------------------
% 
% visualize just one abmn pseudo location and its sensitivity
% 
% ------------------------------------------------------------------------------
nx = 6e2;
nz = 6e2;

x=linspace(0,5,nx);
z=linspace(0,nelectrodes*0.5+1,nz);

[X,Z] = meshgrid(x,z);

dx=x(2)-x(1);
dz=z(2)-z(1);
% ------------------------------------------------------------------------------
figure;
for iexample=1:12
  % choose abmn configuration
  iabmn=randi(nabmn);

  source_p = electrodes(abmn(iabmn,1),:);
  source_n = electrodes(abmn(iabmn,2),:);

  rec_p = electrodes(abmn(iabmn,3),:);
  rec_n = electrodes(abmn(iabmn,4),:);
  % ----------------------------------------------------------------------------
  psi_  = sensitivity_3dDC(X,Z,dx,dz,source_p,source_n,rec_p,rec_n);
  % ----------------------------------------------------------------------------
  subplot(2,6,iexample)
  fancy_imagesc(psi_,x,z)
  colormap(rainbow2_cb(1))
  caxis(5e-3*[-1 1])
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
