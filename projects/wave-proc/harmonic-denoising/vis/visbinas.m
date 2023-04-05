close all
clear
clc
% ------------------------------------------------------------------------------
addpath('../../../pdes/dc-xbore-vis/src/')
% ------------------------------------------------------------------------------
path_read_='../bin/save/';
path_read ='../bin/read/';

path_read_='/media/diegox/7F3D-E672/data/foralles/precis-clu16/round1/rscheckkdensity/save/';
path_read ='/media/diegox/7F3D-E672/data/foralles/precis-clu16/round1/rscheckkdensity/read/';
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
t_=t_+0.01525;
t=(0:(nt-1))*dt;t=t.';
t=t+0.01525;
% ------------------------------------------------------------------------------
%                                    🎼 🎨
% ------------------------------------------------------------------------------
% % because my fft pads the signal with zeros at the end,
% % the fft on heaveside breaks if nt>2^(some power).
% ntpw2=2^nextpow2(nt);
% ntoff = ntpw2-nt;
% dataipsf = [dataips; repmat(dataips(nt,:),[ntoff,1])];
% dataipsf=dataips;
% [dataipsf,f,df] = fourier_rt(dataips,dt);

% ntpw2_=2^nextpow2(nt_);
% ntoff_ = ntpw2_-nt_;
% dataipsf_ = [dataips_; repmat(dataips_(nt_,:),[ntoff_,1])];
% dataipsf_=dataips_;
% [dataipsf_,f_,df_] = fourier_rt(dataips_,dt);
% ------------------------------------------------------------------------------
f_d_f = fft(dataips,[],1);
df = 1/dt/nt;
f = (-nt/2:nt/2-1)*df;
dataipsf = fftshift(f_d_f,1);
dataipsf = dataipsf( ceil(nt/2)+1:nt-1, : );
f = f( ceil(nt/2)+1:nt-1 );

f_d_f = fft(dataips_,[],1);
df = 1/dt/nt_;
f_ = (-nt_/2:nt_/2-1)*df;
dataipsf_ = fftshift(f_d_f,1);
dataipsf_ = dataipsf_( ceil(nt_/2)+1:nt_-1, : );
f_ = f_( ceil(nt_/2)+1:nt_-1 );
% ------------------------------------------------------------------------------
dataipsf = abs(dataipsf) / numel(f);
dataipsf_= abs(dataipsf_) / numel(f_);
% ------------------------------------------------------------------------------
minid=min([min(dataips(:)) min(dataips_(:))]);
maxid=max([max(dataips(:)) max(dataips_(:))]);

minif=min([min(dataipsf(:)) min(dataipsf_(:))]);
maxif=max([max(dataipsf(:)) max(dataipsf_(:))]);
% ------------------------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1],'visible','off');
subplot(1,2,1);
iabmn=1;
loglog(f,dataipsf(:,iabmn),'linewidth',1,'color',rgb(iabmn,:))
hold on;
for iabmn=2:nabmn
loglog(f,dataipsf(:,iabmn),'linewidth',1,'color',rgb(iabmn,:))
end
hold off
axis tight;
axis square;
grid on;
xticks([1,10,100,1000]);
xtickangle(0);
ylim([minif,maxif])
xlabel('Frequency (Hz)')
ylabel('Power (V²/Hz)')
simple_figure();

subplot(1,2,2);
iabmn=1;
loglog(f_,dataipsf_(:,iabmn),'linewidth',1,'color',rgb(iabmn,:))
hold on;
for iabmn=2:nabmn
loglog(f_,dataipsf_(:,iabmn),'linewidth',1,'color',rgb(iabmn,:))
end
hold off
axis tight;
axis square;
grid on;
xticks([1,10,100,1000]);
xtickangle(0);
ylim([minif,maxif])
xlabel('Frequency (Hz)')
ylabel('Power (V²/Hz)')
simple_figure();
% ------------------------------------------------------------------------------
print(gcf,'aspowspec','-dpng','-r350')
% ------------------------------------------------------------------------------
