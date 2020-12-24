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
i_perm = 3000;

b=zeros(24,24);

% 1st column
b(1:8,1:8) = a_perms(i_perm,1);
b(8:16,1:8) = a_perms(i_perm,2);
b(16:24,1:8) = a_perms(i_perm,3);
% 2nd column
b(1:8,8:16) = a_perms(i_perm,4);
b(8:16,8:16) = a_perms(i_perm,5);
b(16:24,8:16) = a_perms(i_perm,6);
% 3rd column
b(1:8,16:24) = a_perms(i_perm,7);
b(8:16,16:24) = a_perms(i_perm,8);
b(16:24,16:24) = a_perms(i_perm,9);

b = b/9;

b_ = image_gaussian_pad(b,kx,ky,'LOW_PASS',nx_pad,ny_pad);
% ------------------------------------------------------------------------------
figure;
subplot(121)
imagesc(b_)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
caxis([1/9 1])
axis square
title('Input')

subplot(122)
imagesc(b)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
caxis([1/9 1])
axis square
title('Output')
% ------------------------------------------------------------------------------
