clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% total number of electrodes
nelectrodes=14; % 32;
% ------------------------------------------------------------------------------
% declare all electrode positions
% ------------------------------------------------------------------------------
% boring setup for borehole stuff
electrodes=[1*ones(nelectrodes/2,1) (linspace(1,nelectrodes*0.5,nelectrodes/2)).'; 4*ones(nelectrodes/2,1) (linspace(1,nelectrodes*0.5,nelectrodes/2)).'];

xmin=min(electrodes(:,1));
xmax=max(electrodes(:,1));
zmin=min(electrodes(:,2));
zmax=max(electrodes(:,2));
% ------------------------------------------------------------------------------
% awesome version for random shit
xmin=1;
xmax=2;
zmin=1;
zmax=2;

electrodes=[xmin+(xmax-xmin)*rand(nelectrodes/2,1)*2 zmin+(zmax-zmin)*rand(nelectrodes/2,1)*2; xmin+(xmax-xmin)*rand(nelectrodes/2,1)*2 zmin+(zmax-zmin)*rand(nelectrodes/2,1)*2];
% ------------------------------------------------------------------------------
% % horizontal borehole
% electrodes=[(linspace(1,nelectrodes*0.5,nelectrodes/2)).' 1*ones(nelectrodes/2,1); (linspace(1,nelectrodes*0.5,nelectrodes/2)).' 4*ones(nelectrodes/2,1) ];
%
% xmin=min(electrodes(:,1));
% xmax=max(electrodes(:,1));
% zmin=min(electrodes(:,2));
% zmax=max(electrodes(:,2));
% ------------------------------------------------------------------------------
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
% repeated    = [2; 4; 6; 10; 11; 15]
% if clusters = [2 4] , [6 10 11 15], then
% klusters    = [2; 4]
% ------------------------------------------------------------------------------
tic;
[klusters_,repeated] = xbore_clusters(abmn,pseud);
toc;
% ------------------------------------------------------------------------------
nrepeat = numel(repeated);
nklu    = size(klusters_,1);
% ------------------------------------------------------------------------------
fprintf('total # of repeated abmn   = %i\n',nrepeat)
fprintf('repeated / total           = %2.2f\n\n',nrepeat / nabmn)
fprintf('total # of clusters        = %i\n\n',nklu)
% ------------------------------------------------------------------------------
% plot electrodes & all pseudo locs (no fancy scheme for repetitions)
pseulocs_plot(electrodes,pseud);
% ------------------------------------------------------------------------------
% plot data in pseudo sections
data2plot = ones(size(abmn,1),1);
units_name='';
cci_path='../../../data/cci_coords/cci';
thics=[40,90,50,0.01];
cmiax=[0,1.5];
pseudats_plot(electrodes,abmn,pseud,klusters_,data2plot,units_name,cci_path,thics,cmiax);
% ------------------------------------------------------------------------------
% random examples of some of the abmn sensitivities
% good for getting an idea of what the survey looks like.
nrows=3;
ncols=6;
prct=3e-4;
sensitiv_plot(electrodes,abmn,pseud,nrows,ncols,prct);
% ------------------------------------------------------------------------------
%
% visualize all mn with one common ab
%
% ------------------------------------------------------------------------------
tic;
s_i_r_d_std = dc_bundle( abmn(:,1:2),ones(nabmn,1),abmn(:,3:4),ones(nabmn,1),zeros(nabmn,1) );
toc;
nsources = size(s_i_r_d_std,2);
% ------------------------------------------------------------------------------
fprintf(' -- just finished bundling abmn -- \n');
fprintf('     there are %i different ab\n\n',nsources);
% ------------------------------------------------------------------------------
nrows=1;
ncols=1;
prct=9e-4;
sensitiv_plots(electrodes,s_i_r_d_std,pseud,nrows,ncols,prct);
% ------------------------------------------------------------------------------
