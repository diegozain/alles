clear
clc
close all
% ..............................................................................
% gigabytes = # of entries * 1e-9 * 8 (for doubles, * 4 for singles)
% 
% for a 1Gb array, we need approx these many entries:
% 
% 1.25e8 ~ 499 * 250501 ~ 500 * 250000 ~ 1e6 * 125 = 1e3 * 1e3 * 125
% 
% so, for a cube of size 1Gb, we have:
% cube_ = zeros(1e3,1e3,125,'double');
% 
% for a plane we could do:
% plane_ = zeros(1e4,12500,'double');
% 
% for a line we do:
% line_ = zeros(1.25e8,1,'double');
% ..............................................................................
nx= 1.25e8;
nx= fix(nx*0.5);
% ..............................................................................
t = linspace(0,2*pi,nx);
dt= t(2)-t(1);
y = sin(t);
% ..............................................................................
% figure;plot(t,y);
% ..............................................................................
% % this takes MORE memory
% yy = struct;
% yy.u = y;

% this takes LESS memory
yy = packer(y);
% ..............................................................................
% clear y;
% ..............................................................................
tic;
integrate_line_(yy,dt);
toc;
% ..............................................................................
% figure;plot(t,yy.u);
% ..............................................................................
delete(yy);
clear
% ..............................................................................