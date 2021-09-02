clc
clear
close all
% ------------------------------------------------------------------------------
addpath('src/');
% ------------------------------------------------------------------------------
na = 2;
nb = 3;
nc = 4;
nd = 5;
% ------------------------------------------------------------------------------
%         _____________
%        /             /|
%       /             / |
%      /_____________/  |
%     |  1           |  |
%     |  2           |  |
% nd  |  3           |  /
%     |  etc         | / nb
%     |______________|/
%           nc
%              .
%              .                  na
%              .
%              .
%         _____________
%        /             /|
%       /             / |
%      /_____________/  |
%     |  1           |  |
%     |  2           |  |
% nd  |  3           |  /
%     |  etc         | / nb
%     |______________|/
%           nc
% ------------------------------------------------------------------------------
% one for loop for four levels
mat4d = zeros(na,nb,nc,nd);
mat4d_= zeros(na*nb*nc*nd,5);
for iabcd = 1:na*nb*nc*nd
  % get d coordinate
  id = mod(iabcd,nd);
  if (id==0)
    id=nd;
  end
  % iabcd = (ia-1)*nb*nc*nd + (ib-1)*nc*nd + (ic-1)*nd + id  ... (*)
  icd = mod(iabcd,nc*nd);
  if (icd==0)
    icd=nc*nd;
  end
  % icd = (ic-1)*nd + id from (*)
  % get c coordinate
  ic = ((icd-id)/nd) + 1;
  % ibcd = (ib-1)*nc*nd + icd
  ibcd = mod(iabcd,nb*nc*nd);
  if (ibcd==0)
    ibcd=nb*nc*nd;
  end
  % get b
  ib = ((ibcd - icd) / (nc*nd)) + 1;
  % get a
  ia = ((iabcd-ibcd)/(nb*nc*nd)) + 1;

  mat4d(ia,ib,ic,id) = iabcd;
  mat4d_(iabcd,:)  = [ia,ib,ic,id,iabcd];
end
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(squeeze(mat4d(1,1,:,:)).')
colormap(rainbow2_cb(1))
caxis([1,na*nb*nc*nd])
xlabel('c')
ylabel('d')
title('slice of 4d array')
simple_figure()

figure;
subplot(1,2,1)
% x y z t
% c b d a
scatter3(mat4d_(1:(nb*nc*nd),3),mat4d_(1:(nb*nc*nd),2),mat4d_(1:(nb*nc*nd),4),80*ones(nb*nc*nd,1),mat4d_(1:(nb*nc*nd),5),'filled')
set(gca,'ZDir','reverse');
axis tight
colormap(rainbow2_cb(1));
caxis([1,na*nb*nc*nd])
colorbar;
xlabel('c')
ylabel('b')
zlabel('d')
title('cube of 4d array')
simple_figure()

subplot(1,2,2)
scatter3(mat4d_((nb*nc*nd+1):(2*nb*nc*nd),3),mat4d_((nb*nc*nd+1):(2*nb*nc*nd),2),mat4d_((nb*nc*nd+1):(2*nb*nc*nd),4),80*ones(nb*nc*nd,1),mat4d_((nb*nc*nd+1):(2*nb*nc*nd),5),'filled')
set(gca,'ZDir','reverse');
axis tight
colormap(rainbow2_cb(1));
caxis([1,na*nb*nc*nd])
colorbar;
xlabel('c')
ylabel('b')
zlabel('d')
title('cube of 4d array')
simple_figure()
% ------------------------------------------------------------------------------
