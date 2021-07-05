clear
close all
clc
% ------------------------------------------------------------------------------
% 
% 1. load data and electrode positions
% 2. bundle all that into cell struct
% 3. compute analytic data one for abmn quadruple
% 
% ------------------------------------------------------------------------------
addpath('src');
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
electr = [posi_dc(:,2) , - posi_dc(:,4)]; % m
nelectr= size(electr,1);
% ------------------------------------------------------------------------------
xmin = min(posi_dc(:,2));
xmax = max(posi_dc(:,2));

ymin = min(-posi_dc(:,4));
ymax = max(-posi_dc(:,4));
% ------------------------------------------------------------------------------
figure;
hold on;
for irecs=1:nelectr
 plot(electr(irecs,1),electr(irecs,2),'k.','markersize',25)
end
% plot(source_p(1),source_p(2),'r.','markersize',25)
% plot(source_n(1),source_n(2),'b.','markersize',25)
hold off;
axis ij
xlim([xmin-2 xmax+2])
ylim([ymin-2 ymax+2])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Electrode positions')
simple_figure();

figure;
scatter(electr(:,1),electr(:,2),50,(1:nelectr),'filled')
colormap(rainbow2_cb(1))
hcb = colorbar;
hcb.Title.String = 'Index #';
axis ij
xlim([xmin-2 xmax+2])
ylim([ymin-2 ymax+2])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Electrode positions')
simple_figure();
% ------------------------------------------------------------------------------
%                  a , b
rec = [data_dc(:,11) , data_dc(:,12)]; % index
src = [data_dc(:,5)  , data_dc(:,8)]; % index

currents_ = data_dc(:,6);

data_dc_o = data_dc(:,13);

std_o = data_dc(:,15);
% ------------------------------------------------------------------------------
% normalize by current magnitude
data_dc_o = data_dc_o ./ currents_;
currents_ = ones(size(currents_));
% ------------------------------------------------------------------------------
% for source-sink #is,
% s_i_r_d_std{is}{1} = [source index , sink index  , current magnitude]
% s_i_r_d_std{is}{2} = [rec + index  , rec - index , obs data , obs std]
s_i_r_d_std = dc_bundle( src,currents_,rec,data_dc_o,std_o );
ns = size(s_i_r_d_std,2);
% ------------------------------------------------------------------------------
fprintf('\n there are %i source-sink pairs\n',ns)
% ------------------------------------------------------------------------------
% choose source-sink index number:
is=1;
% ------------------------------------------------------------------------------
% measuring operator
% of size n_data by 2
% 
% each row is one data-point
% 
% first  column is # of electrode m
% second column is # of electrode n
M = uint32(s_i_r_d_std{is}{2}(:,1:2));
nd= size(M,1);
data_dc = zeros(nd,1);

data_dc_o= s_i_r_d_std{is}{2}(:,3);
source_p = electr(s_i_r_d_std{is}{1}(1),:);
source_n = electr(s_i_r_d_std{is}{1}(2),:);
current_ = s_i_r_d_std{is}{1}(3);
% ------------------------------------------------------------------------------
% 
% two layer solution within the first layer when source is within first layer.
% 
% 
% Electrical Methods & Geophysical Prosource_pecting, Keller. 1927
% pages 107-111.
% https://archive.org/details/electricalmethod00kell/page/106/mode/1up
% ------------------------------------------------------------------------------
u_2l = twolayered_3dDC([sigm0 sigm1 sigm2],h_,electr(:,1),electr(:,2),source_p,source_n,current_);
% ------------------------------------------------------------------------------
for idat = 1:nd
 data_dc(idat) = u_2l(M(idat,1)) - u_2l(M(idat,2));
end
% ------------------------------------------------------------------------------
% 
% 
% 
%                       plot
% 
% 
% 
% ------------------------------------------------------------------------------
nx = 1e2;
nz = 1e2;

x=linspace(min(electr(:,1))-5,max(electr(:,1))+5,nx);
z=linspace(0,max(electr(:,2))+3,nz);
% ------------------------------------------------------------------------------
sigm_=ones(nz,nx);
sigm_(1:binning(z,h_),:) = sigm1;
sigm_((binning(z,h_)+1:nz),:) = sigm2;
% ------------------------------------------------------------------------------
figure;

subplot(1,2,1)
fancy_imagesc(sigm_,x,z)
% colormap(rainbow2_cb(1))
xlabel('Length (m)')
ylabel('Depth (m)')
title('Conductivity (S/m)')
simple_figure()

hold on;
for irecs=1:nelectr
 plot(electr(irecs,1),electr(irecs,2),'k.','markersize',25)
end
plot(source_p(1),source_p(2),'r.','markersize',25)
plot(source_n(1),source_n(2),'b.','markersize',25)
hold off;
% ------------------------------------------------------------------------------
[X,Z] = meshgrid(x,z);
u_2l = twolayered_3dDC([sigm0 sigm1 sigm2],h_,X,Z,source_p,source_n,current_);
u_2l((binning(z,h_)+1:nz),:) = NaN;
% ------------------------------------------------------------------------------
subplot(1,2,2)
fancy_pcolor(u_2l,x,z)
caxis(1e-1*[-max(u_2l(:)) max(u_2l(:))])
colormap(rainbow2_cb(1))
xlabel('Length (m)')
ylabel('Depth (m)')
title('Electric potential (V)')
simple_figure()

hold on;
for irecs=1:nelectr
 plot(electr(irecs,1),electr(irecs,2),'.','markersize',25,'color',[0.4,0.2,0.6])
end
plot(source_p(1),source_p(2),'r.','markersize',25)
plot(source_n(1),source_n(2),'b.','markersize',25)
hold off;
% ------------------------------------------------------------------------------
