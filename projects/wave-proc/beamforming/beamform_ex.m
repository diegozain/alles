close all
clear
clc
% ..............................................................................
% diego domenzain
% fall 2017
% Boise State University
% ..............................................................................
% data bandpass
% ..............................................................................
% below 0.1 Hz is microseismic energy.
% but actually, nothing coherent comes 
% out until 1 Hz.
% above 5 Hz is noise energy.
f_high= 5; % Hz
f_low = 1; % Hz
% ..............................................................................
% beamforming parameters
% ..............................................................................
% individual frequency
fo = 1; % Hz
% frequency range
f_low_ = 1.5; % Hz
f_high_= 3.5; % Hz
fo_ = linspace(f_low_,f_high_,10);
% velocity range
v_min = 2; % Km/s
v_max = 8; % Km/s
v = linspace(v_min,v_max,50);
% angle range
theta = linspace(0,2*pi,100);
% ..............................................................................
path_data = '../../../data/dylan-array/';
name_data = 'array_data_example.mat';
name_locs = 'stationCoordinatesYX_in_km.txt';
% ..............................................................................
receivers = load(strcat(path_data,name_locs));
r = receivers;
r(:,1) = receivers(:,2);
r(:,2) = receivers(:,1);
nr = length(r(:,1));
clear receivers
% ..............................................................................
arr = load(strcat(path_data,name_data));
d = arr.dta;
d = d';
% d is ( nt x nr ) matrix
[nt,nr] = size(d);

% sample rate is 200[Hz]
fs = 200;
fny = fs/2;
dt = 1/fs;
T = (nt-1)*dt;
t = 0:dt:T;
df = fs/nt;
f = (0:nt-1) * df;
% ..............................................................................
figure;
hold on
for i=1:nr
plot(t,d(:,i) + (i-1)*0.2)
end
hold off
axis tight
xlabel('Time (s)')
ylabel('Amplitude')
title('Raw data')
simple_figure();
% ..............................................................................
% detrend and demean
d = detrend(d,1);
d = d-repmat(mean(d,1),nt,1);
% ..............................................................................
% filter
d = filt_gauss(d,dt,f_low,f_high,10);
% ..............................................................................
figure;
subplot(211)
hold on
for i_=1:nr
plot(r(i_,1),r(i_,2),'.','markersize',40)
end
hold off
grid on
axis([-1.2,1.2,-2,0.2])
axis square
xlabel('Length (Km)')
ylabel('Width (Km)')
title('Array of receivers')
simple_figure();

% figure;
subplot(212)
hold on
for i=1:nr
plot(t,d(:,i) + (i-1)*0.2)
end
hold off
axis tight
axis square
xlabel('Time (s)')
set(gca,'yticklabel',[])
title('Filtered data')
simple_figure();
% ..............................................................................
% frequency domain
[d_,f,df] = fourier_rt(d,dt);
% ..............................................................................
% 
%       beamforming velocity vs angle
% 
% ..............................................................................
% individual frequency
b = beamformer_(fo,r,d_,f,v,theta);
% average of an array of frequencies
b = beamformer_thetav(fo_,r,d_,f,v,theta);
% ..............................................................................
b=abs(b).^2;
b=b/max(b(:));
% ..............................................................................
figure;
subplot(122)
fancy_polar(b.',v,theta)
colorbar('off')
colormap(rainbow2(1))
simple_figure();
subplot(121)
imagesc(v,theta,b);
colormap(rainbow2(1))
axis square
xlabel('Velocity (Km/s)')
ylabel('Angle (rad)')
title('Beamform semblance')
simple_figure();
% ..............................................................................
% % for python plotting
% beamformer.b = b;
% beamformer.v = v;
% beamformer.theta = theta;
% beamformer.fo_ = fo_;
% beamformer.f = f;
% beamformer.d_ = d_;
% save(strcat('beamformer.mat'),'beamformer')
% ..............................................................................
% 
%       beamforming slowness vs wavenumber
% 
% ..............................................................................
% build slowness grid
slow_max = 1/v_min;
sx = linspace(-slow_max,slow_max,50);
sy = linspace(-slow_max,slow_max,50);
[sX,sY] = meshgrid(sx,sy);
% ..............................................................................
% individual frequency
b = beamformer(fo,r,d_,f,sX,sY);
% average of frequencies
b = beamformer_sxsy(fo_,r,d_,f,sX,sY);
% ..............................................................................
b=abs(b).^2;
b=b/max(b(:));
% ..............................................................................
figure;
fancy_imagesc(b,sx,sy)
colormap(rainbow2(1))
colorbar('off')
xlabel('Wavenumber in length (1/Km)')
ylabel('Wavenumber in width (1/Km)')
title('Beamform semblance')
simple_figure();
% ..............................................................................
% plot like haney
b = beamformer(fo,r,d_,f,sX,sY);
wo= 2*pi*fo;
beamPlot( b, sX, sY, slow_max, wo );
% ..............................................................................
