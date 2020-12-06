close all
clear
clc
% ..............................................................................
% diego domenzain
% fall 2020
% Colorado School of Mines
% ..............................................................................
to = 3;
fo = 2; % Hz
wo = 2*pi*fo;

co = 0;
bo = 0.8;
% ..............................................................................
% sample rate will be 200 Hz
nt = 3000;
fs = 200;
fny = fs/2;
dt = 1/fs;
T = (nt-1)*dt;
t = 0:dt:T;
t = t.';
% ..............................................................................
wvlet_ = @(to,wo) ( 1-0.5*(wo^2)*(t-to).^2 ) .* exp( -0.25*(wo^2)*(t-to).^2 );
wvlet__= @(to,wo,bo,co) exp(-((t-to).^2)./(bo^2)) .* cos(wo*(t-to) + co);
% ..............................................................................
wvlet = wvlet__(to,wo,bo,co) + wvlet__(2.5*to,2*wo,bo,co) + wvlet__(4*to,4*wo,bo,co) + sin((2*pi/t(nt))*t);
% ..............................................................................
max_=max(wvlet);
min_=min(wvlet);
% ..............................................................................
[stft_, f, gau_] = stft(wvlet, dt, t, 6, 2);
stft_=normali(stft_);
% ..............................................................................
figure('Renderer', 'painters', 'Position', [10 10 600 500]);
subplot(211)
plot(t,wvlet);
ylim([min_-0.1,max_+0.1])
set(gca,'yticklabel',[])
title('Data and STFT')
simple_figure()
subplot(2,1,2)
% fancy_imagesc(abs(stft_),t,f);
fancy_imagesc(abs(stft_(1:binning(f,9),:)),t,f(1:binning(f,9)));
colormap(rainbow2(1))
colorbar('off')
axis xy
axis normal
xlabel('Time (s)')
ylabel('Frequency (Hz)')
simple_figure()
% ..............................................................................
