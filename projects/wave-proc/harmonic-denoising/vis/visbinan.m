close all
clear
clc
% ------------------------------------------------------------------------------
% ambient noise data was recorded with the tx on, but unplugged from the hole.
% however, the (splitter) cable was still connected to the rx.
%
% this means that all data with either m, n, or both m and n in the tx cable were 
% sensitive to ab "injecting" current.
% Moreover, m and n in the tx cable were in contact with the air.
%
% looking at the 50hz denoised data and sorting it by abs(α + β),
% it is clear that the abmn with larger abs(α + β) have
% m, n, or both m and n in the tx cable.
%
% some of these abmn have a clear 18hz signal.
% 
% ------------------------------------------------------------------------------
addpath('../../../pdes/dc-xbore-vis/src/')
% ------------------------------------------------------------------------------
path_read_='../bin/save/';
path_read ='../bin/read/';
% ------------------------------------------------------------------------------
dataips__size= read_bin(strcat(path_read_,'dataips__size'),[3,1],'uint32');
nt_ = dataips__size(1);
nabmn = dataips__size(2);
dataips_ = read_bin(strcat(path_read_,'dataips_'),[nt_*nabmn,1],'single');
dataips_ = reshape(dataips_, [nt_,nabmn]);
bafos_size= read_bin(strcat(path_read_,'bafos_size'),[3,1],'uint32');
nb = bafos_size(1);
nh = bafos_size(2);
nabmn = bafos_size(3);
alphas_ = read_bin(strcat(path_read_,'alphas_'),[nb*nh*nabmn,1],'double');
alphas_ = reshape(alphas_, [nb*nh,nabmn]);
betas_ = read_bin(strcat(path_read_,'betas_'),[nb*nh*nabmn,1],'double');
betas_ = reshape(betas_, [nb*nh,nabmn]);
fos_ = read_bin(strcat(path_read_,'fos_'),[nb*nabmn,1],'double');
fos_ = reshape(fos_, [nb,nabmn]);
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
abmn = read_bin(strcat(path_read,'abmn'),[nabmn,4],'uint32');
dataips_size= read_bin(strcat(path_read,'dataips_size'),[3,1],'uint32');
nt = dataips_size(1);
dataips = read_bin(strcat(path_read,'dataips'),[nt*nabmn,1],'single');
dataips = reshape(dataips, [nt,nabmn]);

dataips_=single(dataips_);
dataips=single(dataips);
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
%                                 NaN 🙅
% ------------------------------------------------------------------------------
ibad=find(isnan(fos_));
fos_(ibad) = [];
alphas_(:,ibad) = [];
betas_(:,ibad) = [];
dataips_(:,ibad) = [];
dataips(:,ibad) = [];
abmn(ibad,:) = [];
nabmn=nabmn-numel(ibad);

ibad=find(isinf(fos_));
fos_(ibad) = [];
alphas_(:,ibad) = [];
betas_(:,ibad) = [];
dataips_(:,ibad) = [];
dataips(:,ibad) = [];
abmn(ibad,:) = [];
nabmn=nabmn-numel(ibad);

ibad=find(fos_< 1e-4);
fos_(ibad) = [];
alphas_(:,ibad) = [];
betas_(:,ibad) = [];
dataips_(:,ibad) = [];
dataips(:,ibad) = [];
abmn(ibad,:) = [];
nabmn=nabmn-numel(ibad);

ibad=find(fos_> 100);
fos_(ibad) = [];
alphas_(:,ibad) = [];
betas_(:,ibad) = [];
dataips_(:,ibad) = [];
dataips(:,ibad) = [];
abmn(ibad,:) = [];
nabmn=nabmn-numel(ibad);
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

% alfabet=alfabet ./ abs(dataips(end,:));

alfabetfos = [alfabet; fos_];
alfabetfos = alfabetfos.';
% ------------------------------------------------------------------------------
[~ , isortab] = sort(alfabetfos(:,1));
alfabetfos = alfabetfos(isortab,:);
fos_=fos_(isortab);
alfabet=alfabet(isortab);
dataips_ = dataips_(:,isortab);
dataips = dataips(:,isortab);
abmn = abmn(isortab,:);
% ------------------------------------------------------------------------------
%
%
%                                     🖌️🖌️🖌️
%
%
% ------------------------------------------------------------------------------
figure;
plot_abmn(abmn)
% ------------------------------------------------------------------------------
rgb=cuatrocolo(nabmn);
rgb=qualitcolor(nabmn);

figure;
subplot(1,3,1)
iabmn=1;
loglog(iabmn,fos_(1),'.','markersize',5,'color',rgb(iabmn,:))
hold on;
for iabmn=2:nabmn
  loglog(iabmn,fos_(iabmn),'.','markersize',5,'color',rgb(iabmn,:))
end
hold off;
axis tight;
axis square;
xticks([1,10,100,1000,10000])
grid on;
ylabel('Frequency (Hz)')
xlabel('# of abmn')
simple_figure()

subplot(1,3,2)
iabmn=1;
loglog(iabmn,alfabet(iabmn),'.','markersize',5,'color',rgb(iabmn,:))
hold on;
for iabmn=2:nabmn
  loglog(iabmn,alfabet(iabmn),'.','markersize',5,'color',rgb(iabmn,:))
end
hold off;
axis tight;
axis square;
xticks([1,10,100,1000,10000])
grid on;
ylabel('α + β')
xlabel('# of abmn')
simple_figure()

subplot(1,3,3)
iabmn=1;
loglog(fos_(iabmn),alfabet(iabmn),'.','markersize',10,'color',rgb(iabmn,:))
hold on;
for iabmn=2:nabmn
  loglog(fos_(iabmn),alfabet(iabmn),'.','markersize',10,'color',rgb(iabmn,:))
end
hold off;
axis tight;
axis square;
grid on;
ylabel('α + β')
xlabel('Frequency (Hz)')
simple_figure()
% ------------------------------------------------------------------------------
%                                  🕐🔌
% ------------------------------------------------------------------------------
dt=2.5e-4;
t_=(0:(nt_-1))*dt;t_=t_.';
t=(0:(nt-1))*dt;t=t.';
% ------------------------------------------------------------------------------
figure;
semilogx(t,dataips(:,1:100),'k')
hold on;
for iabmn=1:100
semilogx(t_,dataips_(:,iabmn),'color',rgb(iabmn,:))
end
hold off;
% ------------------------------------------------------------------------------
mini_=min(dataips_(:));
maxi_=max(dataips_(:));
% ------------------------------------------------------------------------------
%                                    🎼 🎨
% ------------------------------------------------------------------------------
[dataipsf,f,df] = fourier_rt(dataips,dt);
[dataipsf_,f_,df_] = fourier_rt(dataips_,dt);
% ------------------------------------------------------------------------------
figure;
subplot(1,2,1);
loglog(f,abs(dataipsf(:,1:28))/numel(real(dataipsf(:,1:28))),'linewidth',1)%,'color',purpura)
axis tight;
axis square;
grid on;
xticks([1,10,100,1000]);
xtickangle(0);
yticks([1e-8,1e-7,1e-6,1e-5,1e-4])
xlabel('Frequency (Hz)')
ylabel('Power (V²/Hz)')
simple_figure();

subplot(1,2,2);
loglog(f_,abs(dataipsf_(:,1:28))/numel(real(dataipsf_(:,1:28))),'linewidth',1)%,'color',azul)
axis tight;
axis square;
grid on;
xticks([1,10,100,1000]);
xtickangle(0);
yticks([1e-8,1e-7,1e-6,1e-5,1e-4])
xlabel('Frequency (Hz)')
ylabel('Power (V²/Hz)')
simple_figure();
% ------------------------------------------------------------------------------
