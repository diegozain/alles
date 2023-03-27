close all
clear
clc
% ------------------------------------------------------------------------------
addpath('../../../pdes/dc-xbore-vis/src/')
% ------------------------------------------------------------------------------
path_read='../bin/save/';
path_read='E:/data/foralles/precis-clu16/round1/worse/save/';
path_read='E:/data/foralles/precis-clu16/round1/rscheckkdensity/save/';
path_read='E:/data/foralles/noise-clu16/save/';
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
path_read='E:/data/foralles/precis-clu16/round1/worse/read/';
path_read='E:/data/foralles/precis-clu16/round1/rscheckkdensity/read/';
path_read='E:/data/foralles/noise-clu16/read/';
% ------------------------------------------------------------------------------
% abmn_bids = read_bin(strcat(path_read,'abmn_bids'),[nabmn,7],'uint32');

abmn = read_bin(strcat(path_read,'abmn'),[nabmn,4],'uint32');
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
abmn(ibad,:) = [];
nabmn=size(abmn,1);
% ------------------------------------------------------------------------------
%                                      💾
% ------------------------------------------------------------------------------
% save('dataips_','dataips_','-v7.3');
% save('dataips','dataips','-v7.3');
% save('alphas_','alphas_');
% save('betas_','betas_');
% save('fos_','fos_');
% save('abmn','abmn');
% ------------------------------------------------------------------------------
%                                cable 3 was flipped
%
% ------------------------------------------------------------------------------
% 16 10 13 14
u =((32*2+1):(32*3)).';
v =flip(u);
u=(1:32*4).';
v=[(1:32*2).' ; v ; ((32*3+1):32*4).'];
f = translate_u2v(u,v);
for iabmn=1:nabmn
  abmn(iabmn,1) = f(abmn(iabmn,1));
  abmn(iabmn,2) = f(abmn(iabmn,2));
  abmn(iabmn,3) = f(abmn(iabmn,3));
  abmn(iabmn,4) = f(abmn(iabmn,4));
end
% ------------------------------------------------------------------------------
%                                 α + β , 𝐟ₒ
% ------------------------------------------------------------------------------
betas = sum(abs(betas_),1);
alphas= sum(abs(alphas_),1);
alfabet = alphas + betas;
alfabetfos = [alfabet; fos_];
alfabetfos = alfabetfos.';

[~ , isortab] = sort(alfabetfos(:,1));
alfabetfos = alfabetfos(isortab,:);
fos_=fos_(isortab);
alfabet=alfabet(isortab);
dataips_ = dataips_(:,isortab);
dataips = dataips(:,isortab);
abmn = abmn(isortab,:);

alfabetfos(:,1) = (alfabetfos(:,1) - mean(alfabetfos(:,1))) / std(alfabetfos(:,1));
alfabetfos(:,2) = (alfabetfos(:,2) - mean(alfabetfos(:,2))) / std(alfabetfos(:,2));
% ------------------------------------------------------------------------------
%                                     🌂
% ------------------------------------------------------------------------------
% nrow=1; ncol=4;
% nrow=5; ncol=10;
nrow=5; ncol=10;
nclus=nrow*ncol;

% [iclus,clusos] = kmeans(alfabetfos,nclus,'Distance','cosine');
[iclus,clusos] = kmeans(fos_.',nclus,'Distance','sqeuclidean');

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
simple_figure()

subplot(2,2,2)
semilogx(alfabet(find(iclus==1)),'.','markersize',5,'color',rgb(1,:))
hold on;
for iclu=2:nclus
  semilogx(alfabet(find(iclus==iclu)),'.','markersize',5,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
ylabel('|α| + |β|')
xlabel('# of abmn')
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
ylabel('# of abmn')
simple_figure()

subplot(2,2,4)
hold on;
for iclu=1:nclus
  plot(fos_(find(iclus==iclu)),alfabet(find(iclus==iclu)),'.','markersize',10,'color',rgb(iclu,:))
end
hold off;
axis tight;
axis square;
ylabel('|α| + |β|')
xlabel('Frequency (Hz)')
simple_figure()
% ------------------------------------------------------------------------------
figure;
for iclu=1:nclus
  subplot(nrow,ncol,iclu)
  plot_abmn(abmn(find(iclus==iclu),:))
  ylim([1 128])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  axis square;
  xlabel('')
  ylabel('')
end
% ------------------------------------------------------------------------------
figure;
hold on;
for iclu=1:nclus
  plot(iclu,clusos(iclu),'.','markersize',40,'color',rgb(iclu,:))
end
axis tight;
axis square;
grid on;
xlabel('Cluster #')
ylabel('Frequency (Hz)')
simple_figure()
% ------------------------------------------------------------------------------
%                                  🕐🔌
% ------------------------------------------------------------------------------
dt=2.5e-4;
t_=(0:(nt_-1))*dt;t_=t_.';
% t_=t_+0.01525;
t=(0:(nt-1))*dt;t=t.';
% t=t+0.01525;
% ------------------------------------------------------------------------------
ito=binning(t,2);
ito_=binning(t_,2);
t=t(1:ito);
t_=t_(1:ito_);
nt=numel(t);
nt_=numel(t_);
dataips = dataips(1:ito,:);
dataips_ = dataips_(1:ito_,:);
% ------------------------------------------------------------------------------
% dataips = dataips - repmat(dataips(1,:),[nt,1]);
% dataips_= dataips_ - repmat(dataips_(1,:),[nt_,1]);
% ------------------------------------------------------------------------------
mini_=min(dataips_(:));
maxi_=max(dataips_(:));
% ------------------------------------------------------------------------------
figure;
for iclu=1:nclus
  subplot(nrow,ncol,iclu)
  semilogx(t_,dataips_(:,find(iclus==iclu)),'color',rgb(iclu,:))
  ylim([mini_,maxi_])
  % xlim([t(1),2]);
  xlim([dt,2]);
  axis square;
  % xticks([1e-2,1e-1,1])
  % xlabel('Time (sec)')
  % ylabel('Voltage (V)')
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  xlabel('')
  ylabel('')
  simple_figure()
end
% ------------------------------------------------------------------------------
%
%              ▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦
%              ▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦
%              ▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦
%              ▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦
%              ▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦
%
% ------------------------------------------------------------------------------
% addpath('../../../pdes/dc-xbore-vis/src/');
%
% electxyzclu_size= read_bin(strcat(path_read,'electxyzclu_size'),[3,1],'uint32');
% nelect = electxyzclu_size(1);
% electxyzclu = read_bin(strcat(path_read,'electxyzclu'),[nelect,3],'double');
%
% tic;
% pseud = xbore_pseudo(electxyzclu,abmn);
% toc;
%
% klusdata=zeros(nabmn,1);
% for iclu=1:nclus
%   klusdata(find(iclus==iclu)) = iclu;
% end
%
% pseudats_plot3d(electxyzclu,pseud,klusdata,80,'Cluster #');
% colormap(rgb)
%
% pseudats_plot3d(electxyzclu,pseud,alfabet,80,'|α| + |β|');
% rgb_=cytwombly_;
% rgb_=normalizergb(rgb_,min(alfabet),0.02,max(alfabet));
% colormap(rgb_)
% ------------------------------------------------------------------------------
