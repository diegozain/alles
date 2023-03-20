close all
clear
clc
% ------------------------------------------------------------------------------
path_read='../bin/read/';
path_read='E:/data/foralles/precis-clu16/round1/rscheckkdensity/read/';
% path_read='E:/data/foralles/precis-clu16/round1/src9hz/read/';
% path_read='E:/data/foralles/noise-clu16/read/';

dataips_size= read_bin(strcat(path_read,'dataips_size'),[3,1],'uint32');
nt = dataips_size(1);
nabmn = dataips_size(2);
dataips = read_bin(strcat(path_read,'dataips'),[nt,nabmn],'single');
% ------------------------------------------------------------------------------
% abmn = read_bin(strcat(path_read,'abmn'),[nabmn,4],'uint32');
% rschecks = read_bin(strcat(path_read,'rschecks'),[nabmn,1],'uint32');
% ikeep=find(rschecks);
% dataips = dataips(:,ikeep);
% abmn = abmn(ikeep,:);
% nabmn = numel(ikeep);
% dataips_size(2) = nabmn;
% save_bin(strcat(path_read,'abmn'),abmn,'uint32');
% save_bin(strcat(path_read,'dataips_size'),dataips_size,'uint32');
% save_bin(strcat(path_read,'dataips'),dataips,'single');
% ------------------------------------------------------------------------------
path_read='../bin/save/';
path_read='E:/data/foralles/precis-clu16/round1/rscheckkdensity/save/';
% path_read='E:/data/foralles/precis-clu16/round1/src9hz/save/';
% path_read='E:/data/foralles/noise-clu16/save/';

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
dt=2.5e-4;
t=(0:(nt-1))*dt;
t=t.';
t=t+0.01525;
t_=t(1:nt_);
% ------------------------------------------------------------------------------
% ito=binning(t,2);
% ito_=binning(t_,2);
% t=t(1:ito);
% t_=t_(1:ito_);
% nt=numel(t);
% nt_=numel(t_);
% dataips = dataips(1:ito,:);
% dataips_ = dataips_(1:ito_,:);
% ------------------------------------------------------------------------------
%
%                                🎵🎵🎵🎵
%
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%
%                               🎨🎨🎨🎨
%
% ------------------------------------------------------------------------------
rgb=cuatrocolo(nabmn);
% ------------------------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 0.7],'visible','off');
subplot(1,3,1)
hold on;
for iabmn=1:nabmn
  plot(iabmn,fos_(iabmn),'.','markersize',25,'colo',rgb(iabmn,:));
end
hold off;
axis square;
axis tight;
xlabel('abmn #')
ylabel('fo (Hz)')
simple_figure()

subplot(1,3,2)
hold on;
for iabmn=1:nabmn
  plot(alphas_(:,iabmn),'.-','markersize',25,'colo',rgb(iabmn,:));
end
hold off;
axis square;
axis tight;
xlabel('Harmonic #')
ylabel('α')
simple_figure()

subplot(1,3,3)
hold on;
for iabmn=1:nabmn
plot(betas_(:,iabmn),'.-','markersize',25,'colo',rgb(iabmn,:));
end
hold off;
axis square;
axis tight;
xlabel('Harmonic #')
ylabel('β')
simple_figure()
% ------------------------------------------------------------------------------
print(gcf,'noise9hz-1','-dpng','-r350')
% ------------------------------------------------------------------------------
ifig=1;
figure('units','normalized','outerposition',[0 0 0.6 1],'visible','off');
for iabmn_=1:2
  iabmn=randi(nabmn);

  subplot(2,2,ifig)
  semilogx(t,1e3*dataips(:,iabmn),'linewidth',2)
  hold on;
  semilogx(t_,1e3*dataips_(:,iabmn),'linewidth',2)
  hold off;
  axis square;
  axis tight;
  grid on;
  xlim([dt,2])
  xticks([1e-3,1e-2,1e-1,1])
  % xlim([0,2])
  % xticks([0,1,2])
  % xlim([0,4])
  % xticks([1e-3,1e-2,1e-1,1])
  xlabel('Time (sec)')
  ylabel('Voltage (mV)')
  % ylabel('Ampers (mA)')
  simple_figure()

  ifig=ifig+1;

  subplot(2,2,ifig)
  plot(t,1e3*dataips(:,iabmn),'linewidth',2)
  hold on;
  plot(t_,1e3*dataips_(:,iabmn),'linewidth',2)
  hold off;
  axis square;
  axis tight;
  grid on;
  xlim([1,2])
  xticks([1,1.5,2])
  % xlim([1,4])
  % xticks([1,2,4])
  xlabel('Time (sec)')
  ylabel('Voltage (mV)')
  % ylabel('Ampers (mA)')
  simple_figure()

  ifig=ifig+1;
end
% ------------------------------------------------------------------------------
print(gcf,'noise9hz-2','-dpng','-r350')
% ------------------------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 0.7 0.7],'visible','off');
subplot(1,2,1)
semilogx(t,1e3*dataips(:,1),'linewidth',1,'color',rgb(1,:))
hold on;
for iabmn=2:nabmn
semilogx(t,1e3*dataips(:,iabmn),'linewidth',1,'color',rgb(iabmn,:))
end
hold off;
axis square;
axis tight;
grid on;
colormap(rgb);
xlim([t(1),2])
xticks([1e-1,1])
% xlim([0,4])
% xticks([1e-3,1e-2,1e-1,1])
xlabel('Time (sec)')
ylabel('Voltage (mV)')
% ylabel('Ampers (mA)')
simple_figure()

subplot(1,2,2)
semilogx(t_,1e3*dataips_(:,1),'linewidth',1,'color',rgb(1,:))
hold on;
for iabmn=2:nabmn
semilogx(t_,1e3*dataips_(:,iabmn),'linewidth',1,'color',rgb(iabmn,:))
end
hold off;
axis square;
axis tight;
grid on;
colormap(rgb);
xlim([t_(1),2])
xticks([1e-1,1])
% xlim([0,4])
% xticks([1e-3,1e-2,1e-1,1])
xlabel('Time (sec)')
ylabel('Voltage (mV)')
% ylabel('Ampers (mA)')
simple_figure()
% ------------------------------------------------------------------------------
print(gcf,'noise9hz-3','-dpng','-r350')
% ------------------------------------------------------------------------------
[dataipspw,f,df] = fourier_rt(dataips,dt);
dataipspw = abs(dataipspw) / numel(f);

[dataipspw_,f_,df] = fourier_rt(dataips_,dt);
dataipspw_ = abs(dataipspw_) / numel(f_);

mini=min([min(dataipspw(:)) min(dataipspw_(:))]);
maxi=max([max(dataipspw(:)) max(dataipspw_(:))]);

figure('units','normalized','outerposition',[0 0 0.7 0.7],'visible','off');
subplot(1,2,1)
loglog(f,dataipspw(:,1),'linewidth',1,'color',rgb(1,:))
hold on;
for iabmn=2:nabmn
loglog(f,dataipspw(:,iabmn),'linewidth',1,'color',rgb(iabmn,:))
end
hold off;
axis square;
axis tight;
grid on;
xticks([1,10,100,1000]);
ylim([mini,maxi])
xlabel('Frequency (Hz)')
ylabel('Power (V²/Hz)')
% ylabel('Power (A²/Hz)')
simple_figure()

subplot(1,2,2)
loglog(f_,dataipspw_(:,1),'linewidth',1,'color',rgb(1,:))
hold on;
for iabmn=2:nabmn
loglog(f_,dataipspw_(:,iabmn),'linewidth',1,'color',rgb(iabmn,:))
end
hold off;
axis square;
axis tight;
grid on;
xticks([1,10,100,1000]);
ylim([mini,maxi])
xlabel('Frequency (Hz)')
ylabel('Power (V²/Hz)')
% ylabel('Power (A²/Hz)')
simple_figure()
% ------------------------------------------------------------------------------
print(gcf,'noise9hz-4','-dpng','-r350')
% ------------------------------------------------------------------------------
dataips = dataips - repmat(dataips(1,:),[nt,1]);
dataips_= dataips_ - repmat(dataips_(1,:),[nt_,1]);

mini=min([min(dataips(:)) min(dataips_(:))]);
maxi=max([max(dataips(:)) max(dataips_(:))]);
% ------------------------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 0.7 0.7],'visible','off');
subplot(1,2,1)
semilogx(t,1e3*dataips(:,1),'linewidth',1,'color',rgb(1,:))
hold on;
for iabmn=2:nabmn
semilogx(t,1e3*dataips(:,iabmn),'linewidth',1,'color',rgb(iabmn,:))
end
hold off;
axis square;
axis tight;
grid on;
colormap(rgb);
xlim([t(1),2])
xticks([1e-1,1])
% xlim([0,4])
% xticks([1e-3,1e-2,1e-1,1])
ylim([mini,maxi]*1e3)
xlabel('Time (sec)')
ylabel('Voltage anchored at 0 (mV)')
% ylabel('Amperage anchored at 0 (mA)')
simple_figure()

subplot(1,2,2)
semilogx(t_,1e3*dataips_(:,1),'linewidth',1,'color',rgb(1,:))
hold on;
for iabmn=2:nabmn
semilogx(t_,1e3*dataips_(:,iabmn),'linewidth',1,'color',rgb(iabmn,:))
end
hold off;
axis square;
axis tight;
grid on;
colormap(rgb);
xlim([t(1),2])
xticks([1e-1,1])
% xlim([0,4])
% xticks([1e-3,1e-2,1e-1,1])
ylim([mini,maxi]*1e3)
xlabel('Time (sec)')
ylabel('Voltage anchored at 0 (mV)')
% ylabel('Amperage anchored at 0 (mA)')
simple_figure()
% ------------------------------------------------------------------------------
print(gcf,'noise9hz-5','-dpng','-r350')
% ------------------------------------------------------------------------------
close all;
clear;
% ------------------------------------------------------------------------------
