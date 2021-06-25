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
% 1. load data
% ------------------------------------------------------------------------------
% posi_dc=load('riprap/fwd_posi.rip');
% data_dc=load('riprap/fwd_data.rip');

posi_dc=load('riprap/fwd_posi.rap');
data_dc=load('riprap/fwd_data.rap');
% ------------------------------------------------------------------------------
% halfspace parameters
h_=22; % m
sigm0=0;
sigm1=1;
sigm2=1e2;
% ------------------------------------------------------------------------------
%  electrodes in meters
% it is implicitly assumed that y is constant for all electrode positions
electrodes = [posi_dc(:,2) , - posi_dc(:,4)]; % m
nelectrodes= size(electrodes,1);
% ------------------------------------------------------------------------------
xmin = min(posi_dc(:,2));
xmax = max(posi_dc(:,2));

zmin = min(-posi_dc(:,4));
zmax = max(-posi_dc(:,4));
% ------------------------------------------------------------------------------
figure;
scatter(electrodes(:,1),electrodes(:,2),50,(1:nelectrodes),'filled')
colormap(rainbow2_cb(1))
hcb = colorbar;
hcb.Title.String = 'Index #';
axis ij
xlim([xmin-2 xmax+2])
ylim([zmin-2 zmax+2])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Electrode positions')
simple_figure();
% ------------------------------------------------------------------------------
rec = [data_dc(:,11) , data_dc(:,12)]; % index m & n
src = [data_dc(:,5)  , data_dc(:,8)]; % index a & b

currents_ = data_dc(:,6);
data_dc_o = data_dc(:,13);
std_o     = data_dc(:,15);

clear data_dc;
% ------------------------------------------------------------------------------
% normalize by current magnitude
data_dc_o = data_dc_o ./ currents_;
currents_ = ones(size(currents_));
% ------------------------------------------------------------------------------
% 2. bundle
% ------------------------------------------------------------------------------
% for source-sink #is,
% s_i_r_d_std{is}{1} = [source index , sink index  , current magnitude]
% s_i_r_d_std{is}{2} = [rec + index  , rec - index , obs data , obs std]
s_i_r_d_std = dc_bundle( src,currents_,rec,data_dc_o,std_o );
ns = size(s_i_r_d_std,2);

clear data_dc_o;
% ------------------------------------------------------------------------------
fprintf('\n there are %i source-sink pairs\n',ns)
% ------------------------------------------------------------------------------
% 3. compute analytic data 
% ------------------------------------------------------------------------------
% choose source-sink index number:
data_dc  =[];
data_dc_o=[];
for is=1:ns
  % ----------------------------------------------------------------------------
  % measuring operator
  % of size n_data by 2
  % 
  % each row is one data-point
  % 
  % first  column is # of electrode m
  % second column is # of electrode n
  M = uint32(s_i_r_d_std{is}{2}(:,1:2));
  nd= size(M,1);
  data_dc_ = zeros(nd,1);

  data_dc_o_= s_i_r_d_std{is}{2}(:,3);
  source_p  = electrodes(s_i_r_d_std{is}{1}(1),:);
  source_n  = electrodes(s_i_r_d_std{is}{1}(2),:);
  current_  = s_i_r_d_std{is}{1}(3);
  % ----------------------------------------------------------------------------
  % 
  % two layer solution within the first layer when source is within first layer.
  % 
  % 
  % Electrical Methods & Geophysical Prosource_pecting, Keller. 1927
  % pages 107-111.
  % https://archive.org/details/electricalmethod00kell/page/106/mode/1up
  % ----------------------------------------------------------------------------
  u_2l = twolayered_3dDC([sigm0 sigm1 sigm2],h_,electrodes(:,1),electrodes(:,2),source_p,source_n,current_);
  % ----------------------------------------------------------------------------
  for idat = 1:nd
    data_dc_(idat) = u_2l(M(idat,1)) - u_2l(M(idat,2));
  end
  data_dc   = [data_dc; data_dc_];
  data_dc_o = [data_dc_o; data_dc_o_];
end
% ------------------------------------------------------------------------------
error_ = sqrt(sum((data_dc - data_dc_o).^2)) / sqrt(sum((data_dc_o).^2));
% ------------------------------------------------------------------------------
fprintf('\n      the error is %f \n\n',error_)
% ------------------------------------------------------------------------------
% 
% 
% 
%                       plot
% 
% 
% 
% ------------------------------------------------------------------------------
figure;
hold on;
plot(data_dc_o,'k.','markersize',40);
plot(data_dc,'r.','markersize',25);
hold off;
legend({'observed','synthetic'})
xlabel('Index #')
ylabel('Data (V)')
simple_figure()
% ------------------------------------------------------------------------------
% 
% 
% 
%                       pseudo-sections part
% 
% 
% 
% ------------------------------------------------------------------------------
abmn = [src , rec];
nabmn=size(abmn,1);
% ------------------------------------------------------------------------------
% get all pseudo locations
tic;
pseud= xbore_pseudo(electrodes,abmn);
toc;
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
figure;
subplot(1,2,1)
hold on;
for ie=1:nelectrodes
plot(electrodes(ie,1),electrodes(ie,2),'k.','markersize',40)
end
for iabmn=1:nabmn
plot(pseud(iabmn,1),pseud(iabmn,2),'r.','markersize',20)
% plot(pseud(iabmn,1),pseud(iabmn,2),'.')
end
hold off;
axis ij
xlim([xmin-2 xmax+2])
ylim([zmin-2 zmax+2])
xlabel('Length (m)')
ylabel('Depth (m)')
title('All AB.MN pseudo locations')
simple_figure()

subplot(1,2,2)
hold on;
for ie=1:nelectrodes
plot(electrodes(ie,1),electrodes(ie,2),'.','markersize',40,'color',[0.5,0.5,0.5])
end
axis ij
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
xlim([xmin-2 xmax+2])
ylim([zmin-2 zmax+2])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Repeated locations')
simple_figure()
% ------------------------------------------------------------------------------
% 
% visualize just one abmn pseudo location and its sensitivity
% 
% ------------------------------------------------------------------------------
nx = 6e2;
nz = 6e2;

x=linspace(xmin-2,xmax+2,nx);
z=linspace(zmin-2,zmax+2,nz);

[X,Z] = meshgrid(x,z);

dx=x(2)-x(1);
dz=z(2)-z(1);
% ------------------------------------------------------------------------------
figure;
for iexample=1:4
  % choose abmn configuration
  iabmn=randi(nabmn);

  source_p = electrodes(abmn(iabmn,1),:);
  source_n = electrodes(abmn(iabmn,2),:);

  rec_p = electrodes(abmn(iabmn,3),:);
  rec_n = electrodes(abmn(iabmn,4),:);
  % ----------------------------------------------------------------------------
  psi_  = sensitivity_3dDC(X,Z,dx,dz,source_p,source_n,rec_p,rec_n);
  % ----------------------------------------------------------------------------
  subplot(2,2,iexample)
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