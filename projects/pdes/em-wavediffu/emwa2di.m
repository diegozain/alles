close all
clear
clc
% ------------------------------------------------------------------------------
fprintf('\n\n\n                                   🤔🤔🤔 \n            is the propagation diffusive 🌼 and the media conductive ⚡? \n                                   🤔🤔🤔 \n\n\n\n')
% ------------------------------------------------------------------------------
% is EM wave propagation asymptotically equal to diffusion?
% if this happens then yes,
%
%                            √(t² - R²/v²) >> 2ε/σ
%                        ⟹             R² << v² ⋅ (t² - 4ε²/σ²)
%
% t is time that has passed and,
% R is distance from source to receiver.
%
% "Diffusion of electromagnetic fields into a two-dimensional earth: a finite difference approach"
% Oristaglio and Hohmann, 1984.
% ------------------------------------------------------------------------------
c = 3e8; % m/s
epso = 8.8e-12; % F/m
% ------------------------------------------------------------------------------
epsr = 50;
eps = epsr * epso;
v = c/sqrt(epso);
% ------------------------------------------------------------------------------
% t = logspace(-11,-9,1e3); % s
% t = logspace(-11,-8,1e3); % s
t = logspace(-9,-3,1e3); % s
sig = logspace(-3,0,1e3); % S/m

[sig_,t_] = meshgrid(sig,t);
% ------------------------------------------------------------------------------
%                             R² << v² ⋅ (t² - 4ε²/σ²)
%
%                         let R_ = √(v² ⋅ (t² - 4ε²/σ²))
% ------------------------------------------------------------------------------
% if R_ is large and positive, then it is diffusion
% if R_ is small and positive, then it is on the edge of being a wave
% if R_ is negative, then it is a wave 💯
% 🚨🚨 large & small compared to closest offset rx.
% ------------------------------------------------------------------------------
R_ = v^2 * (t_.^2 - (4*eps^2)./(sig_.^2));
R_(find(R_(:)<0)) = NaN;
R_ = sqrt(R_);
% ------------------------------------------------------------------------------
%                                  🎨🎨🎨🎨
% ------------------------------------------------------------------------------
figure;
contour(log10(sig),log10(t),log10(R_),50,'linewidth',10);
axis square;
axis ij;
colormap(rainbow2_cb(1));
colorbar;
ylabel('Time (lg s)')
xlabel('Conductivity (lg S/m)')
title('√(v²t² - 4v²ε²/σ²) (lg m)')
simple_figure()
% ------------------------------------------------------------------------------
% is the media a poor conductor or a good conductor?
%
% σ << ωε poor conductor
% σ >> ωε good conductor
%
%     or,
%
% ε >> σ/ω poor conductor
% ε << σ/ω good conductor
% ------------------------------------------------------------------------------
sig = logspace(-3,0,1e3); % S/m
freqs = logspace(-10,9,1e3);
omegs = 2*pi*freqs;

[sig_,omegs_] = meshgrid(sig,omegs);
% ------------------------------------------------------------------------------
% if epsis is large, then the media is a good conductor
epsis = sig_ ./ omegs_;
epsis = epsis / epso;
epsis(find(epsis(:)<1)) = NaN;
% ------------------------------------------------------------------------------
%                                  🎨🎨🎨🎨
% ------------------------------------------------------------------------------
figure;
contour(log10(sig),log10(omegs),log10(epsis),50,'linewidth',10);
axis square;
axis ij;
colormap(rainbow2_cb(1));
colorbar;
ylabel('Frequency (lg rad)')
xlabel('Conductivity (lg S/m)')
title('σ/ω (lg - )')
simple_figure()
% ------------------------------------------------------------------------------
