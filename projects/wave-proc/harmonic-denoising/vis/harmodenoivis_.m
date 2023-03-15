close all
clear
clc
% ------------------------------------------------------------------------------
addpath('../../../pdes/dc-xbore-vis/src/')
% ------------------------------------------------------------------------------
path_read='../bin/save/';

dataips__size= read_bin(strcat(path_read,'dataips__size'),[3,1],'uint32');
nt_ = dataips__size(1);
nabmn = dataips__size(2);
dataips_ = read_bin(strcat(path_read,'dataips_'),[nt_*nabmn,1],'single');
dataips_ = reshape(dataips_, [nt_,nabmn]);
bafos_size= read_bin(strcat(path_read,'bafos_size'),[3,1],'uint32');
nb = bafos_size(1);
nh = bafos_size(2);
nabmn = bafos_size(3);
alphas_ = read_bin(strcat(path_read,'alphas_'),[nb*nh*nabmn,1],'double');
alphas_ = reshape(alphas_, [nb*nh,nabmn]);
betas_ = read_bin(strcat(path_read,'betas_'),[nb*nh*nabmn,1],'double');
betas_ = reshape(betas_, [nb*nh,nabmn]);
fos_ = read_bin(strcat(path_read,'fos_'),[nb*nabmn,1],'double');
fos_ = reshape(fos_, [nb,nabmn]);
% ------------------------------------------------------------------------------
path_read='../bin/read/';

abmn = read_bin(strcat(path_read,'abmn'),[nabmn,4],'uint32');
dataips_size= read_bin(strcat(path_read,'dataips_size'),[3,1],'uint32');
nt = dataips_size(1);
dataips = read_bin(strcat(path_read,'dataips'),[nt*nabmn,1],'single');
dataips = reshape(dataips, [nt,nabmn]);

dataips_=single(dataips_);
dataips=single(dataips);
% ------------------------------------------------------------------------------
ibad=find(isnan(fos_));
fos_(ibad) = [];
alphas_(:,ibad) = [];
betas_(:,ibad) = [];
dataips_(:,ibad) = [];
dataips(:,ibad) = [];
abmn(ibad,:) = [];
nabmn=size(abmn,1);
% ------------------------------------------------------------------------------
% save('dataips_','dataips_','-v7.3');
% save('dataips','dataips','-v7.3');
% save('alphas_','alphas_');
% save('betas_','betas_');
% save('fos_','fos_');
% save('abmn','abmn');
% ------------------------------------------------------------------------------
betas = sum(abs(betas_),1);
alphas= sum(abs(alphas_),1);
% ------------------------------------------------------------------------------
nabm=zeros(nabmn,4);
nabm(:,1)=abmn(:,4);
nabm(:,2)=abmn(:,1);
nabm(:,3)=abmn(:,2);
nabm(:,4)=abmn(:,3);

[ab,mn,ab_mn,iabmn_sort] = abmn_bundle(nabm);
nabm=nabm(iabmn_sort,:);
fos_=fos_(iabmn_sort);
alphas=alphas(iabmn_sort);
betas=betas(iabmn_sort);
dataips_=dataips_(:,iabmn_sort);
dataips=dataips(:,iabmn_sort);
% ------------------------------------------------------------------------------
nclus=12;
[iclus,cluso] = kmeans(fos_.',nclus,'Distance','sqeuclidean');

clussizes = zeros(nclus,1);
for iclu=1:nclus
  clussizes(iclu) = numel(find(iclus==iclu));
end
% ------------------------------------------------------------------------------
rgb=cuatrocolo(nclus);
% ------------------------------------------------------------------------------
figure;
subplot(3,3,4)
semilogx(fos_(find(iclus==1)),'.','markersize',5,'color',rgb(1,:))
hold on;
for iclu=2:nclus
semilogx(fos_(find(iclus==iclu)),'.','markersize',5,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
ylabel('Frequency (Hz)')
xlabel('# of abmn')

subplot(3,3,5)
semilogy(1,clussizes(1),'.','markersize',20,'color',rgb(1,:))
hold on;
for iclu=2:nclus
semilogy(iclu,clussizes(iclu),'.','markersize',20,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
xlabel('Cluster #')
ylabel('# of abmn')

subplot(3,3,6)
semilogy(fos_(find(iclus==1)),clussizes(1),'.','markersize',20,'color',rgb(1,:))
hold on;
for iclu=2:nclus
  semilogy(fos_(find(iclus==iclu)),clussizes(iclu),'.','markersize',20,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
xlabel('Frequency (Hz)')
ylabel('# of abmn')
% ------------------------------------------------------------------------------
subplot(3,3,1)
semilogx(alphas(find(iclus==1)),'.','markersize',5,'color',rgb(1,:))
hold on;
for iclu=2:nclus
semilogx(alphas(find(iclus==iclu)),'.','markersize',5,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
ylabel('α')
xlabel('# of abmn')

subplot(3,3,3)
semilogy(alphas(find(iclus==1)),clussizes(1),'.','markersize',20,'color',rgb(1,:))
hold on;
for iclu=2:nclus
  semilogy(alphas(find(iclus==iclu)),clussizes(iclu),'.','markersize',20,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
xlabel('α')
ylabel('# of abmn')
% ------------------------------------------------------------------------------
subplot(3,3,7)
semilogx(betas(find(iclus==1)),'.','markersize',5,'color',rgb(1,:))
hold on;
for iclu=2:nclus
semilogx(betas(find(iclus==iclu)),'.','markersize',5,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
ylabel('β')
xlabel('# of abmn')

subplot(3,3,9)
semilogy(betas(find(iclus==1)),clussizes(1),'.','markersize',20,'color',rgb(1,:))
hold on;
for iclu=2:nclus
  semilogy(betas(find(iclus==iclu)),clussizes(iclu),'.','markersize',20,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
xlabel('β')
ylabel('# of abmn')
% ------------------------------------------------------------------------------
figure;
for iclu=1:nclus
  % subplot(5,10,iclu)
  subplot(3,4,iclu)
  plot_abmn(nabm(find(iclus==iclu),:))
  % legend({'n','a','b','m'})
  ylim([1 128])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  axis square;
  xlabel('')
  ylabel('')
end
% ------------------------------------------------------------------------------
dt=2.5e-4;
t_=(0:(nt_-1))*dt;t_=t_.';t_=t_+0.01525;
t=(0:(nt-1))*dt;t=t.';t=t+0.01525;
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
