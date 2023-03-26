close all
clear
clc
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
% ------------------------------------------------------------------------------
% [~ , isortab] = sort(alfabetfos(:,1));
% alfabetfos = alfabetfos(isortab,:);
% fos_=fos_(isortab);
% alfabet=alfabet(isortab);
% dataips_ = dataips_(:,isortab);
% dataips = dataips(:,isortab);
% ------------------------------------------------------------------------------
%
%
%                                     🖌️🖌️🖌️
%
%
% ------------------------------------------------------------------------------
rgb=cuatrocolo(nabmn);

% ------------------------------------------------------------------------------
%                                  🕐🔌
% ------------------------------------------------------------------------------
dt=2.5e-4;
t_=(0:(nt_-1))*dt;t_=t_.';
t=(0:(nt-1))*dt;t=t.';
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
ylabel('Power (A²/Hz)')
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
ylabel('Power (A²/Hz)')
simple_figure();
% ------------------------------------------------------------------------------
