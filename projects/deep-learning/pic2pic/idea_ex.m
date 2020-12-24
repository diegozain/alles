clear
close all
clc
% ------------------------------------------------------------------------------
%
%
%         this script is just an example of what should happen.
%
%
% ------------------------------------------------------------------------------
a=(1:9).';
a_perms = perms(a);
% ------------------------------------------------------------------------------
% set smoothing
% larger number = less smoothing
kx=0.1;
ky=kx;
nx_pad = 36;
ny_pad = 36;
% ------------------------------------------------------------------------------
% choose permutation
iperm = 3000;

b=zeros(24,24);

b(1:8,1:8) = a_perms(iperm,1);
b(8:16,1:8) = a_perms(iperm,2);
b(16:24,1:8) = a_perms(iperm,3);

b(1:8,8:16) = a_perms(iperm,4);
b(8:16,8:16) = a_perms(iperm,5);
b(16:24,8:16) = a_perms(iperm,6);

b(1:8,16:24) = a_perms(iperm,7);
b(8:16,16:24) = a_perms(iperm,8);
b(16:24,16:24) = a_perms(iperm,9);

b_ = image_gaussian_pad(b,kx,ky,'LOW_PASS',nx_pad,ny_pad);
% ------------------------------------------------------------------------------
figure;
subplot(121)
imagesc(b_)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
caxis([1 9])
axis square
title('Input')

subplot(122)
imagesc(b)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
caxis([1 9])
axis square
title('Output')
% ------------------------------------------------------------------------------
