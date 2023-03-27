close all
clear
clc
% ------------------------------------------------------------------------------
addpath('../../../pdes/dc-xbore-vis/src/')
% ------------------------------------------------------------------------------
path_read='../bin/save/';
% path_read='E:/data/foralles/klus4clara/save/';
% path_read='E:/data/foralles/precis-clu16/round1/src9hz/save/';
% ------------------------------------------------------------------------------
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
% path_read='E:/data/foralles/klus4clara/read/';
% path_read='E:/data/foralles/precis-clu16/round1/src9hz/read/';
% ------------------------------------------------------------------------------
dataips_size= read_bin(strcat(path_read,'dataips_size'),[3,1],'uint32');
nt = dataips_size(1);
dataips = read_bin(strcat(path_read,'dataips'),[nt*nabmn,1],'single');
dataips = reshape(dataips, [nt,nabmn]);

dataips_=single(dataips_);
dataips=single(dataips);
% ------------------------------------------------------------------------------
%                                 NaN 🙅
% ------------------------------------------------------------------------------
ibad=find(isnan(fos_));
fos_(ibad) = [];
alphas_(:,ibad) = [];
betas_(:,ibad) = [];
dataips_(:,ibad) = [];
dataips(:,ibad) = [];
% ------------------------------------------------------------------------------
%                                      💾
% ------------------------------------------------------------------------------
% save('dataips_','dataips_','-v7.3');
% save('dataips','dataips','-v7.3');
% save('alphas_','alphas_');
% save('betas_','betas_');
% save('fos_','fos_');
% ------------------------------------------------------------------------------
%                                  α + β , 𝐟ₒ
% ------------------------------------------------------------------------------
betas = sum(abs(betas_),1);
alphas= sum(abs(alphas_),1);
alfabet = alphas + betas;

alfabet=alfabet ./ abs(dataips(end,:));

alfabetfos = [alfabet; fos_];
alfabetfos = alfabetfos.';

[~ , isortab] = sort(alfabetfos(:,1));
alfabetfos = alfabetfos(isortab,:);
fos_=fos_(isortab);
alfabet=alfabet(isortab);
dataips_ = dataips_(:,isortab);
dataips = dataips(:,isortab);

alfabetfos(:,1) = (alfabetfos(:,1) - mean(alfabetfos(:,1))) / std(alfabetfos(:,1));
alfabetfos(:,2) = (alfabetfos(:,2) - mean(alfabetfos(:,2))) / std(alfabetfos(:,2));
% ------------------------------------------------------------------------------
%                                     🌂
% ------------------------------------------------------------------------------
% nrow=1; ncol=5;
% nrow=5; ncol=10;
nrow=1; ncol=3;
nclus=nrow*ncol;

[iclus,clusos] = kmeans(alfabetfos,nclus,'Distance','cosine');
% [iclus,clusos] = kmeans(fos_.',nclus,'Distance','sqeuclidean');

clussizes = zeros(nclus,1);
for iclu=1:nclus
  clussizes(iclu) = numel(find(iclus==iclu));
end
% ------------------------------------------------------------------------------
%
%
%                                     🖌️🖌️🖌️
%
%
% ------------------------------------------------------------------------------
rgb=cuatrocolo(nclus);

figure;
subplot(2,2,1)
semilogx(fos_(find(iclus==1)),'.','markersize',15,'color',rgb(1,:))
hold on;
for iclu=2:nclus
  semilogx(fos_(find(iclus==iclu)),'.','markersize',15,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
ylabel('Frequency (Hz)')
xlabel('# of ab')
simple_figure()

subplot(2,2,2)
semilogx(alfabet(find(iclus==1)),'.','markersize',15,'color',rgb(1,:))
hold on;
for iclu=2:nclus
  semilogx(alfabet(find(iclus==iclu)),'.','markersize',15,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
ylabel('|α| + |β|')
xlabel('# of ab')
simple_figure()

subplot(2,2,3)
semilogy(1,clussizes(1),'.','markersize',40,'color',rgb(1,:))
hold on;
for iclu=2:nclus
  semilogy(iclu,clussizes(iclu),'.','markersize',40,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
xlabel('Cluster #')
ylabel('# of ab')
simple_figure()

subplot(2,2,4)
hold on;
for iclu=1:nclus
  plot(fos_(find(iclus==iclu)),alfabet(find(iclus==iclu)),'.','markersize',15,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
ylabel('|α| + |β|')
xlabel('Frequency (Hz)')
simple_figure()
% ------------------------------------------------------------------------------
dt=2.5e-4;
t_=(0:(nt_-1))*dt;t_=t_.';
t=(0:(nt-1))*dt;t=t.';
% ------------------------------------------------------------------------------
mini_=min(dataips_(:));
maxi_=max(dataips_(:));
% ------------------------------------------------------------------------------
figure;
for iclu=1:nclus
  subplot(nrow,ncol,iclu)
  semilogx(t_,1e3*dataips_(:,find(iclus==iclu)),'color',rgb(iclu,:))
  ylim([mini_,maxi_]*1e3)
  xticks([1e-2,1e-1,1])
  axis square;
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  xlabel('')
  ylabel('')
  simple_figure()
end
% ------------------------------------------------------------------------------
%                                    🎼 🎨
% ------------------------------------------------------------------------------
[dataipsf,f,df] = fourier_rt(dataips,dt);
[dataipsf_,f_,df_] = fourier_rt(dataips_,dt);
% ------------------------------------------------------------------------------
figure;
subplot(1,2,1);
loglog(f,abs(dataipsf)/numel(real(dataipsf)),'linewidth',1)%,'color',purpura)
axis tight;
axis square;
grid on;
xticks([1,10,100,1000]);
xtickangle(0);
yticks([1e-8,1e-7,1e-6,1e-5,1e-4])
xlabel('Frequency (Hz)')
ylabel('Power (A²/Hz)')
simple_figure();

subplot(1,2,2);
loglog(f_,abs(dataipsf_)/numel(real(dataipsf_)),'linewidth',1)%,'color',azul)
axis tight;
axis square;
grid on;
xticks([1,10,100,1000]);
xtickangle(0);
yticks([1e-8,1e-7,1e-6,1e-5,1e-4])
xlabel('Frequency (Hz)')
ylabel('Power (A²/Hz)')
simple_figure();
% ------------------------------------------------------------------------------
