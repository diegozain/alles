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
a=(1:4).';
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
i_perm = 15;

b=zeros(24,24);

% 1st column
b(1:12,1:12) = a_perms(i_perm,1);
b(13:24,1:12) = a_perms(i_perm,2);
% 2nd column
b(1:12,13:24) = a_perms(i_perm,3);
b(13:24,13:24) = a_perms(i_perm,4);

b = b/4;

b_ = image_gaussian_pad(b,kx,ky,'LOW_PASS',nx_pad,ny_pad);
% ------------------------------------------------------------------------------
figure;
subplot(121)
imagesc(b_)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
caxis([1/4 1])
axis square
title('Input')

subplot(122)
imagesc(b)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
caxis([1/4 1])
axis square
title('Output')
% ------------------------------------------------------------------------------