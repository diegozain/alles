close all
clear
clc
% ..............................................................................
% diego domenzain
% fall 2017
% Boise State University
% ..............................................................................
prompt = '\n   do you want an average or mono frequency beamforming? (a or m): ';
ave_mono = input(prompt,'s');
% ..............................................................................
% synthetic parameters
% ..............................................................................
vo = 4; % Km/s
fo = 2; % Hz
wo = 2*pi*fo;
% ..............................................................................
% beamforming parameters
% ..............................................................................
% frequency range
f_low_ = 0.5; % Hz
f_high_= 5; % Hz
fo_ = linspace(f_low_,f_high_,100);
% velocity range
v_min = 2.5; % Km/s
v_max = 8; % Km/s
v = linspace(v_min,v_max,50);
% angle range
theta = linspace(0,2*pi,100);
% ..............................................................................
path_data = '../../../data/dylan-array/';
name_locs = 'stationCoordinatesYX_in_km.txt';
% ..............................................................................
receivers = load(strcat(path_data,name_locs));
r = receivers;
r(:,1) = receivers(:,2);
r(:,2) = receivers(:,1);
nr = length(r(:,1));
clear receivers
% ..............................................................................
sources = zeros(nr,2);

% sources(:,1) = -70;
% sources(:,2) = r(:,2);

sources(:,2) = -70;
sources(:,1) = r(:,1);
% ..............................................................................
% sample rate will be 200 Hz
nt = 7200;
fs = 200;
fny = fs/2;
dt = 1/fs;
T = (nt-1)*dt;
t = 0:dt:T;
t = t.';
% ..............................................................................
% synth data
% ..............................................................................
% wavelet
wvlet = @(to) ( 1-0.5*(wo^2)*(t-to).^2 ) .* exp( -0.25*(wo^2)*(t-to).^2 );
d = zeros(nt,nr);

for ir = 1:nr
 to_ = norm(sources(ir,:)-r(ir,:))/vo;
 d(:,ir) = wvlet(to_);
end
% ..............................................................................
figure('Renderer', 'painters', 'Position', [10 10 400 700]);
subplot(211)
hold on
for ir=1:nr
 plot(r(ir,1),r(ir,2),'.','markersize',40)
% plot(sources(ir,1),sources(ir,2),'rp','markersize',10,'MarkerFaceColor','red')
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
plot(t,d(:,i) + (i-1)*2)
end
hold off
axis tight
axis square
xlabel('Time (s)')
set(gca,'yticklabel',[])
title('Data')
simple_figure();
% ..............................................................................
% frequency domain
[d_,f,df] = fourier_rt(d,dt);
% ..............................................................................
% 
%       beamforming velocity vs angle
% 
% ..............................................................................
if strcmp(ave_mono,'m')
 % individual frequency
 b = beamformer_(fo,r,d_,f,v,theta);
elseif strcmp(ave_mono,'a')
 % % average of an array of frequencies
 b = beamformer_thetav(fo_,r,d_,f,v,theta);
end
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
if strcmp(ave_mono,'m')
 % individual frequency
 b = beamformer(fo,r,d_,f,sX,sY);
elseif strcmp(ave_mono,'a')
 % average of frequencies
 b = beamformer_sxsy(fo_,r,d_,f,sX,sY);
end
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
beamPlot( b, sX, sY, slow_max, wo );
% ..............................................................................