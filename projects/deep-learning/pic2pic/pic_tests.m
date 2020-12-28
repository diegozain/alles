clear
close all
clc
% ------------------------------------------------------------------------------
% this script produces pictures to test, different from the training set.
% 
% the output pictures are just an example of what should happen.
%
%
% ------------------------------------------------------------------------------
% set smoothing
% larger number = less smoothing
kx=0.1;
ky=kx;
nx_pad = 36;
ny_pad = 36;
% ------------------------------------------------------------------------------
figure;
% ------------------------------------------------------------------------------
% a square in the middle
c1 = ones(24,24);
c1(7:18,7:18) = 4;

c1 = c1/4;

subplot(151)
imagesc(c1)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
caxis([1/4 1])
axis square
title('Square')
% ------------------------------------------------------------------------------
% a diagonal matrix
c2 = ones(24,24);
c2 = triu(c2);
c2(c2==0) = 4;

c2 = c2/4;

c3 = image_gaussian_pad(c2,kx,ky,'LOW_PASS',nx_pad,ny_pad);

subplot(152)
imagesc(c2)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
caxis([1/4 1])
axis square
title('Diagonal')

subplot(153)
imagesc(c3)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
caxis([1/4 1])
axis square
title('Smooth diagonal')
% ------------------------------------------------------------------------------
% a circle in the middle
kx = 6;
ky = 6;

X = 0:24-1;
Y = 0:24-1;
x = 12;
y = 12;
[XY,YX] = meshgrid(X,Y);
c4 = ((XY-x)/kx).^2 + ((YX-y)/ky).^2;
c4 = exp(-c4);

c4 = (c4+1)*0.5;

subplot(154)
imagesc(c4)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
caxis([1/4 1])
axis square
title('Smooth circle')
% ------------------------------------------------------------------------------
% random but smoothed
c5 = randi(4,24,24);

c5 = c5/4;

c5 = image_gaussian_pad(c5,2e-2*kx,2e-2*ky,'LOW_PASS',nx_pad,ny_pad);

subplot(155)
imagesc(c5)
simple_figure()
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap(rainbow2(1))
% caxis([1/4 1])
axis square
title('Random smooth')
% ------------------------------------------------------------------------------
c =zeros(24,24,5);

c(:,:,1) = c1;
c(:,:,2) = c2;
c(:,:,3) = c3;
c(:,:,4) = c4;
c(:,:,5) = c5;
% ------------------------------------------------------------------------------
% save
save('../../../data/pic2pic/c','c')
% ------------------------------------------------------------------------------
