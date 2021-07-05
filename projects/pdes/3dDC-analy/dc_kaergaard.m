clear
close all
clc
% ------------------------------------------------------------------------------
% 
% 1. load data and electrode positions
% 2. bundle all that into cell struct
% 3. compute analytic data & error for all abmn quadruples
%  
% 4. get pseudo locations and plot that too
% 
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% 1. load positions
% ------------------------------------------------------------------------------
posi_dc=load('riprap/kaergaard_posi.rip');
% ------------------------------------------------------------------------------
%  electrodes in meters: x,y,z
electrodes = [posi_dc(:,2) , posi_dc(:,3) , - posi_dc(:,5)]; % m
nelectrodes= size(electrodes,1);
% ------------------------------------------------------------------------------
xmin = min(posi_dc(:,2));
xmax = max(posi_dc(:,2));

ymin = min(posi_dc(:,3));
ymax = max(posi_dc(:,3));

zmin = min(-posi_dc(:,5));
zmax = max(-posi_dc(:,5));

clear posi_dc;
% ------------------------------------------------------------------------------
% the dimensions of these x & y are enormoous, so I'm gonna make them human 
electrodes(:,1) = electrodes(:,1)-xmin;
electrodes(:,2) = electrodes(:,2)-ymin;

xmax = xmax-xmin;
ymax = ymax-ymin;
xmin = 0;
ymin = 0;
% ------------------------------------------------------------------------------
xyunique = unique(electrodes(:,1:2),'rows');
% ------------------------------------------------------------------------------
figure;
% 3d plot
subplot(1,2,1)
scatter3(electrodes(:,1),electrodes(:,2),electrodes(:,3),50,(1:nelectrodes),'filled')
colormap(rainbow2_cb(1))
hcb = colorbar;
hcb.Title.String = 'Index #';
set(gca, 'YDir','reverse')
set(gca, 'ZDir','reverse')
xlim([xmin-2 xmax+2])
ylim([ymin-2 ymax+2])
xlabel('Length (m)')
ylabel('Width (m)')
zlabel('Depth (m)')
title('Electrode positions xyz')
simple_figure();

% 2d-xy plot
subplot(1,2,2)
plot(xyunique(:,1),xyunique(:,2),'k.','markersize',30)
axis ij
xlim([xmin-2 xmax+2])
ylim([ymin-2 ymax+2])
xlabel('Length (m)')
ylabel('Width (m)')
title('Electrode positions xy')
simple_figure();
% ------------------------------------------------------------------------------
% 2. load data
% ------------------------------------------------------------------------------
data_dc=load('riprap/kaergaard_data.rip');
% ------------------------------------------------------------------------------
rec = [data_dc(:,11) , data_dc(:,12)]; % index m & n
src = [data_dc(:,5)  , data_dc(:,8)]; % index a & b

currents_ = data_dc(:,6);
data_dc_o = data_dc(:,14); % apparent resistivity. NOT voltage
std_o     = data_dc(:,15);

clear data_dc;
% ------------------------------------------------------------------------------
% this data is so noisy
irhoa_neg = find(data_dc_o<=0);

rec(irhoa_neg,:) = [];
src(irhoa_neg,:) = [];
currents_(irhoa_neg) = [];
data_dc_o(irhoa_neg) = [];
std_o(irhoa_neg) = [];
% ------------------------------------------------------------------------------
fprintf('\n I just removed %i bad data points. bye-bye.\n',numel(irhoa_neg))
% ------------------------------------------------------------------------------
% 3. bundle
% ------------------------------------------------------------------------------
% % for source-sink #is,
% % s_i_r_d_std{is}{1} = [source index , sink index  , current magnitude]
% % s_i_r_d_std{is}{2} = [rec + index  , rec - index , obs data , obs std]
% s_i_r_d_std = dc_bundle( src,currents_,rec,data_dc_o,std_o );
% ns = size(s_i_r_d_std,2);
% fprintf('\n there are %i source-sink pairs,',ns)
% ------------------------------------------------------------------------------
abmn = [src , rec];
abmn = uint32(abmn);
nabmn= size(abmn,1);
% ------------------------------------------------------------------------------
clear irhoa_neg src currents_ rec std_o
% ------------------------------------------------------------------------------
fprintf('\n there are %i abmn quadruples\n\n',nabmn)
% ------------------------------------------------------------------------------
% 4. choose boreholes
% ------------------------------------------------------------------------------
% indexes of electrodes per borehole
% I got these manually, couldnt bother.
% b1 = 1:30;
% b2 = 31:58;
% b3 = 59:90;
% b4 = 91:115;
% b5 = 116:144;
% b6 = 145:168;
% b7 = 169:192;
% b8 = 193:216;
% b9 = 217:240;
% ------------------------------------------------------------------------------
% boreholes 3 & 5 look pretty centered. Borehole 5 is the middle one.
b3 = 59:90;
b5 = 116:144;
ib = [b3 , b5];
electrodes__=electrodes(ib,:);
electrodes__(:,2)=[];
nelectrodes__=size(electrodes__,1);

electrodes_xz =electrodes;
electrodes_xz(:,2)=[];
% ------------------------------------------------------------------------------
% get abmn
nabmn_=0;
for iabmn=1:nabmn
  a=abmn(iabmn,1);
  b=abmn(iabmn,2);
  m=abmn(iabmn,3);
  n=abmn(iabmn,4);
  if (a> 58 && a< 91 && b> 58 && b< 91 && m> 115 && m< 145 && n> 115 && n< 145)
    nabmn_=nabmn_+1;
  end
end

abmn_=zeros(nabmn_,4);
data_dc=zeros(nabmn_,1);
iabmn_=0;
for iabmn=1:nabmn
  a=abmn(iabmn,1);
  b=abmn(iabmn,2);
  m=abmn(iabmn,3);
  n=abmn(iabmn,4);
  if (a> 58 && a< 91 && b> 58 && b< 91 && m> 115 && m< 145 && n> 115 && n< 145)
    iabmn_=iabmn_+1;
    abmn_(iabmn_,:) = abmn(iabmn,:);
    data_dc(iabmn_) = data_dc_o(iabmn);
  end
end
% ------------------------------------------------------------------------------
% get all pseudo locations
tic;
pseud= xbore_pseudo(electrodes_xz,abmn_);
toc;
% ------------------------------------------------------------------------------
% get clusters of repeated elements
tic;
[klusters_,repeated] = xbore_clusters(abmn_,pseud);
toc;
% ------------------------------------------------------------------------------
nrepeat = numel(repeated);
nklu = size(klusters_,1);
% ------------------------------------------------------------------------------
fprintf('\ntotal # of repeated abmn   = %i\n',nrepeat)
fprintf('repeated / total           = %2.2f\n\n',nrepeat / nabmn_)
fprintf('total # of clusters        = %i\n\n',nklu)
% ------------------------------------------------------------------------------
figure;
hold on;
for ie=1:nelectrodes__
plot(electrodes__(ie,1),electrodes__(ie,2),'k.','markersize',40)
end
for iabmn=1:nabmn_
plot(pseud(iabmn,1),pseud(iabmn,2),'r.','markersize',20)
% plot(pseud(iabmn,1),pseud(iabmn,2),'.')
end
hold off;
axis ij;
axis image;
xlim([1 5])
ylim([zmin-1 zmax+1])
xlabel('Length (m)')
ylabel('Depth (m)')
title('All AB.MN pseudo locations')
simple_figure()
% ------------------------------------------------------------------------------
% 
% 
% visualize the data
% 
% 
% ------------------------------------------------------------------------------
figure;
hold on;
% plot electrodes
scatter(electrodes__(:,1),electrodes__(:,2),150,'k','filled')
% plot all pseudo positions (including repetitions)
scatter(pseud(:,1),pseud(:,2),160,log10(data_dc),'filled')
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
  
  pseud_ = pseud(klusters_{iklu},:) + 0.2*dense_circle;
  data_dc_=data_dc(klusters_{iklu});
  
  scatter(pseud_(:,1),pseud_(:,2),20,log10(data_dc_),'filled')
  % ----------------------------------------------------------------------------
end
hold off;
colormap(rainbow2_cb(1))
hcb = colorbar;
hcb.Title.String = 'log10(Ohm.m)';
axis ij;
axis image;
xlim([1 5])
ylim([zmin-1 zmax+1])
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

x=linspace(1,5,nx);
z=linspace(0,zmax+2,nz);

[X,Z] = meshgrid(x,z);

dx=x(2)-x(1);
dz=z(2)-z(1);
% ------------------------------------------------------------------------------
figure;
for iexample=1:10
  % choose abmn configuration
  iabmn=randi(nabmn_);

  source_p = electrodes_xz(abmn_(iabmn,1),:);
  source_n = electrodes_xz(abmn_(iabmn,2),:);

  rec_p = electrodes_xz(abmn_(iabmn,3),:);
  rec_n = electrodes_xz(abmn_(iabmn,4),:);
  % ----------------------------------------------------------------------------
  psi_  = sensitivity_3dDC(X,Z,dx,dz,source_p,source_n,rec_p,rec_n);
  % ----------------------------------------------------------------------------
  subplot(2,5,iexample)
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
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
