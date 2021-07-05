clc
clear
close all
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% total number of electrodes
nelectrodes=32; % 32; 
% declare all electrode positions
electrodes=[1*ones(nelectrodes/2,1) , zeros(nelectrodes/2,1), 9+(linspace(1,nelectrodes*0.6,nelectrodes/2)).'; 4*ones(nelectrodes/2,1) , zeros(nelectrodes/2,1), 9+(linspace(1,nelectrodes*0.6,nelectrodes/2)).'];
% # of electrodes that will be Tx
nTx=nelectrodes/2;
% ------------------------------------------------------------------------------
xmin = min(electrodes(:,1));
zmin = min(electrodes(:,3));

xmax = max(electrodes(:,1));
zmax = max(electrodes(:,3));
% ------------------------------------------------------------------------------
% get all abmn pairs
abmn = xbore_getall(nelectrodes,nTx);
nabmn= size(abmn,1);
% ------------------------------------------------------------------------------
figure;
plot(electrodes(:,1),electrodes(:,3),'k.','markersize',40);
axis ij
xlim([xmin-1, xmax+1])
ylim([0, zmax+1])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Survey electrodes')
simple_figure()
% ------------------------------------------------------------------------------
% the idea is to run the "analytical" solution using the 1d forward model 
% in AarhusInv.
% we compare this "analytical" solution to the 3d forward model living in the 
% Line inversion.
% ------------------------------------------------------------------------------
load('riprap/forward_1d.fwr')
load('riprap/forward_2layer.rap')
% ------------------------------------------------------------------------------
data_1d=forward_1d(:,13);
data_2l=forward_2layer(:,14);
% ------------------------------------------------------------------------------
i1dnans=find(isnan(data_1d));
abmn_1dnans = abmn(i1dnans,:);

abmn(i1dnans,:)=[];
data_1d(i1dnans)=[];
data_2l(i1dnans)=[];

figure;hold on;plot(data_1d,'r.','markersize',30);plot(data_2l,'k.','markersize',25);
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
