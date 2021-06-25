clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
%          3d DC analytic solution 
% ------------------------------------------------------------------------------
%      x   z
s_p = [10 ; 0];
s_m = [20 ; 0];

current_ = 1; % A

sigm=1; % S/m

% h =10; % m
h_=10;% m
sigm0=0;
sigm1=1;
sigm2=100;
n_accuracy = 50;
% ------------------------------------------------------------------------------
% 
%                       analytic potential everywhere
% 
% ------------------------------------------------------------------------------
nx = 1e2;
nz = 1e2;

x=linspace(0,30,nx);
z=linspace(0,30,nz);

[X,Z] = meshgrid(x,z);
% ------------------------------------------------------------------------------
s_p = [x(binning(x,s_p(1))) ; z(binning(z,s_p(2)))];
s_m = [x(binning(x,s_m(1))) ; z(binning(z,s_m(2)))];

s_pim = [x(binning(x,s_p(1))) ; -z(binning(z,s_p(2)))];
s_mim = [x(binning(x,s_m(1))) ; -z(binning(z,s_m(2)))];
% ------------------------------------------------------------------------------
dist_p = sqrt((X - s_p(1)).^2 + (Z - s_p(2)).^2);
dist_m = sqrt((X - s_m(1)).^2 + (Z - s_m(2)).^2);

dist_pim = sqrt((X - s_pim(1)).^2 + (Z - s_pim(2)).^2);
dist_mim = sqrt((X - s_mim(1)).^2 + (Z - s_mim(2)).^2);
% ------------------------------------------------------------------------------
u_ho = current_/(4*pi*sigm) * ( 1./dist_p - 1./dist_m + 1./dist_pim - 1./dist_mim  );
% ------------------------------------------------------------------------------
subplot(3,4,[1,2,5,6])
fancy_imagesc(u_ho,x,z)
colormap(rainbow2_cb(1))
xlabel('Length (m)')
ylabel('Depth (m)')
title('Homogeneous (V)')
simple_figure()
% ------------------------------------------------------------------------------
subplot(3,4,[9,10])
semilogy(x,abs(u_ho(1,:)),'k.-','markersize',20)
grid on
ylim([5e-3 1])
xlabel('Length (m)')
ylabel('| Surface potential | (V)')
simple_figure()
% ------------------------------------------------------------------------------
% 
% two layer solution within the first layer when source is within first layer.
% 
% 
% Electrical Methods & Geophysical Prosource_pecting, Keller. 1927
% pages 107-111.
% https://archive.org/details/electricalmethod00kell/page/106/mode/1up
% ------------------------------------------------------------------------------
u_2l = twolayered_3dDC([sigm0 sigm1 sigm2],h_,X,Z,s_p,s_m,current_);
% ------------------------------------------------------------------------------
subplot(3,4,[3,4,7,8])
fancy_imagesc(u_2l,x,z)
colormap(rainbow2_cb(1))
xlabel('Length (m)')
ylabel('Depth (m)')
title('Two-layered (V)')
simple_figure()
% ------------------------------------------------------------------------------
subplot(3,4,[11,12])
semilogy(x,abs(u_2l(1,:)),'k.-','markersize',20)
grid on
ylim([5e-3 1])
xlabel('Length (m)')
% ylabel('| Surface potential | (V)')
simple_figure()
% ------------------------------------------------------------------------------
figure;
semilogy(x,abs(u_2l(1,:)),'k.-','markersize',20)
grid on
ylim([5e-3 1])
xlabel('Length (m)')
ylabel('| Surface potential | (V)')
simple_figure()
% ------------------------------------------------------------------------------
