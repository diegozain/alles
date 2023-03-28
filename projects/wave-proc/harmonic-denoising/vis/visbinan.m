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
% path_read_='../bin/save/';
% path_read ='../bin/read/';

% path_read_='E:data/foralles/noise-clu16/9hz/save/';
% path_read ='E:data/foralles/noise-clu16/9hz/read/';

% path_read_='E:data/foralles/noise-clu16/50hz/save/';
% path_read ='E:data/foralles/noise-clu16/50hz/read/';

path_read_='D:data/an50hz/save/';
path_read ='D:data/an50hz/read/';
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

abmn_bids = read_bin(strcat(path_read,'abmn_bids'),[nabmn,7],'uint32');

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

% ------------------------------------------------------------------------------
iabcab=[];
ibacabls=[];
iokcabls=[];
for iabmn=1:nabmn
  abid=abmn_bids(iabmn,5);
  mid=abmn_bids(iabmn,6);
  nid=abmn_bids(iabmn,7);
  if (abid==mid && abid==nid)
    iabcab = [iabcab; iabmn];
  end
  if (abid==mid && abid~=nid)
    ibacabls = [ibacabls; iabmn];
  end
  if (abid==nid && abid~=mid)
    ibacabls = [ibacabls; iabmn];
  end
  if (abid~=mid && abid~=nid)
    iokcabls = [iokcabls; iabmn];
  end
end

fos_ = fos_(ibacab);
alphas_ = alphas_(:,ibacab);
betas_ = betas_(:,ibacab);
dataips_ = dataips_(:,ibacab);
dataips = dataips(:,ibacab);
abmn = abmn(ibacab,:);
nabmn=numel(ibacab);
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

ibad=find(fos_< 0.5);
fos_(ibad) = [];
alphas_(:,ibad) = [];
betas_(:,ibad) = [];
dataips_(:,ibad) = [];
dataips(:,ibad) = [];
abmn(ibad,:) = [];
nabmn=nabmn-numel(ibad);

ibad=find(fos_> 2e3-1);
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
% ------------------------------------------------------------------------------
% alfabet=alfabet ./ abs(dataips(end,:));

% ibad=find(alfabet> 0.1);
% fos_(ibad) = [];
% alphas_(:,ibad) = [];
% betas_(:,ibad) = [];
% alfabet(ibad) = [];
% dataips_(:,ibad) = [];
% dataips(:,ibad) = [];
% abmn(ibad,:) = [];
% nabmn=nabmn-numel(ibad);
% ------------------------------------------------------------------------------
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
plot_abmn(abmn);
print(gcf,'anabmn','-dpng','-r350')

mini=min(alfabet);
maxi=max(alfabet);
% ------------------------------------------------------------------------------
rgb=cuatrocolo(nabmn);
rgb=qualitcolor(nabmn);

figure('units','normalized','outerposition',[0 0 1 0.7]);
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
% yticks([0.1,1,10])
yticks([0.1,10,100])
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
xticks([1,10,100,1000,10000]);
ylim([mini,maxi]);
grid on;
ylabel('|α| + |β|')
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
% xticks([0.1,1,10]);
xticks([0.1,10,100])
ylim([mini,maxi]);
grid on;
ylabel('|α| + |β|')
xlabel('Frequency (Hz)')
simple_figure()
% ------------------------------------------------------------------------------
print(gcf,'analfabetfos','-dpng','-r350')
% ------------------------------------------------------------------------------
%
%                                  🕐🔌
%
% ------------------------------------------------------------------------------
dt=2.5e-4;
t_=(0:(nt_-1))*dt;t_=t_.';
t=(0:(nt-1))*dt;t=t.';
ito=nt; ito_=nt_;
% ------------------------------------------------------------------------------
%                                   🔪✂️
% ------------------------------------------------------------------------------
% ito=binning(t,2);
% ito_=binning(t_,2);
% dataips=dataips(1:ito,:);
% dataips_=dataips_(1:ito_,:);
% ------------------------------------------------------------------------------
%                                    🎼 🎨
% ------------------------------------------------------------------------------
[dataipsf,f,df] = fourier_rt(dataips,dt);
[dataipsf_,f_,df_] = fourier_rt(dataips_,dt);

dataipsf = abs(dataipsf) / numel(f);
dataipsf_= abs(dataipsf_) / numel(f_);
% ------------------------------------------------------------------------------
minid=min([min(dataips(:)) min(dataips_(:))]);
maxid=max([max(dataips(:)) max(dataips_(:))]);

minif=min([min(dataipsf(:)) min(dataipsf_(:))]);
maxif=max([max(dataipsf(:)) max(dataipsf_(:))]);
% %{
% ------------------------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1],'visible','off');
subplot(1,2,1)
iabmn=1;
semilogx(t(1:ito),dataips(:,iabmn),'color',rgb(iabmn,:))
hold on;
for iabmn=2:nabmn
semilogx(t(1:ito),dataips(:,iabmn),'color',rgb(iabmn,:))
end
hold off;
axis tight;
axis square;
grid on;
xticks([1e-4,1e-3,1e-2,1e-1,1]);
xtickangle(0);
ylim([minid,maxid])
xlabel('Time (sec)')
ylabel('Voltage (V)')
simple_figure();

subplot(1,2,2)
iabmn=1;
semilogx(t_(1:ito_),dataips_(:,iabmn),'color',rgb(iabmn,:))
hold on;
for iabmn=2:nabmn
semilogx(t_(1:ito_),dataips_(:,iabmn),'color',rgb(iabmn,:))
end
hold off;
axis tight;
axis square;
grid on;
xticks([1e-4,1e-3,1e-2,1e-1,1]);
xtickangle(0);
ylim([minid,maxid])
xlabel('Time (sec)')
ylabel('Voltage (V)')
simple_figure();
% ------------------------------------------------------------------------------
print(gcf,'antimedat','-dpng','-r350')
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
print(gcf,'anpowspec','-dpng','-r350')
% ------------------------------------------------------------------------------
%}
