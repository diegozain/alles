clc
clear
close all
% ------------------------------------------------------------------------------
load('ip_data')
dt=ip_data.dt;
d=ip_data.d;

d(1:3)=[];

nt=numel(d);
t=(0:dt:(dt*(nt-1))).';
% ------------------------------------------------------------------------------
% d = window_mean(d,20);
% ------------------------------------------------------------------------------
% % tukey cannot be activated in the filter!!
% [d,~]    = filt_gauss(d,dt,-3e1,3e1,8);
% ------------------------------------------------------------------------------
[d_,f,df] = fourier_rt(d,dt);
% ------------------------------------------------------------------------------
figure;

subplot(2,1,1)
hold on
plot(t,d,'k-','linewidth',3);
plot(t,zeros(nt,1),'--','color',[0.5,0.5,0.5],'linewidth',2)
hold off
xlabel('Time (s)')
ylabel('Voltage (V)')
simple_figure()

subplot(2,1,2)
plot(f,normali(abs(d_)),'k-','linewidth',3);
xlabel('Frequency (Hz)')
ylabel('Normalized power')
simple_figure()
% ------------------------------------------------------------------------------
figure;

y = differentiate_line(normali(abs(d_)),df);

hold on
plot(f,y,'k-','linewidth',3);
plot(f,zeros(numel(f),1),'--','color',[0.5,0.5,0.5],'linewidth',2)
hold off
xlabel('Frequency (Hz)')
ylabel('Derivative of power')
simple_figure()
% ------------------------------------------------------------------------------
