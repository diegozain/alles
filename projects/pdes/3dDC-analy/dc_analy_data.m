clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% halfspace parameters
h_=7; % m
sigm0=0;
sigm1=1;
sigm2=1e2;
% ------------------------------------------------------------------------------
% sources
%      x   z
source_p = [5 , 5]; % m
source_n = [5 , 6]; % m
% current
current_ = 1; % A
%  electrodes in meters
electr = [5 1; 5 2; 5 3; 5 4; 5 5; 5 6; 15 1; 15 2; 15 3; 15 4; 15 5; 15 6]; % m
nelectr= size(electr,1);
% ------------------------------------------------------------------------------
% measuring operator
% of size n_data by 2
% 
% each row is one data-point
% 
% first  column is # of electrode m
% second column is # of electrode n
M = [3 4; 5 6; 7 8; 9 10; 11 12];
M = uint32(M);
nd= size(M,1);
data_dc = zeros(nd,1);
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
caxis([-0.5 0.5])
colormap(rainbow2_cb(1))
xlabel('Length (m)')
ylabel('Depth (m)')
title('Electric potential (V)')
simple_figure()

hold on;
for irecs=1:nelectr
 plot(electr(irecs,1),electr(irecs,2),'k.','markersize',25)
end
plot(source_p(1),source_p(2),'r.','markersize',25)
plot(source_n(1),source_n(2),'b.','markersize',25)
hold off;
% ------------------------------------------------------------------------------
