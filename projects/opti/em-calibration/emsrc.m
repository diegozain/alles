close all
clear
clc
% ------------------------------------------------------------------------------
%
%             450
%          .------.
%          |      |
%          |      |
%  .-------.      .---------.
%  0      10⁵   5.5×10⁵     10⁶   ⟵  ns
%  0      100    550       1000   ⟵  μs
%  0     1e-4   5.5e-4     1e-3   ⟵  s
%
% ------------------------------------------------------------------------------
% but then i was like 🤔 maybe we could do ricker no?
% so the ricker ones are within the step -- 5 × 10³ Hz = 5 kHz
% ------------------------------------------------------------------------------
%                                   👷👷👷👷👷⌚
% ------------------------------------------------------------------------------
% total time
ttot = 1e-3; % s
% 1s = 10⁶ μs
dt = 1e-6 * 10; % 10; % s
nt = fix(ttot / dt);
t = ((0:(nt-1))*dt).';
% center time
tce = 325;
% μs ⟶ s
tce = tce*1e-6; % s
% ------------------------------------------------------------------------------
%                                   🍘 Ricker 🍘
% ------------------------------------------------------------------------------
% from "Frequencies of the Ricker wavelet". Yanghua Wang.
fo = 5e3; % Hz
wo = 2*pi*fo;
tau = 1 / pi / fo;
src = ( 1-0.5*(wo^2)*(t-tce).^2 ) .* exp( -0.25*(wo^2)*(t-tce).^2 );
% ------------------------------------------------------------------------------
%                              ∂t ( 🍘 Ricker 🍘 )
% ------------------------------------------------------------------------------
% fo = 5e3; % Hz
% tau = 1 / pi / fo;
% src = -(t-tce).*exp( -( ((t-tce).^2) ./ (tau^2) ) );
% ------------------------------------------------------------------------------
%                              🚶 step function 🚶
% ------------------------------------------------------------------------------
% % width
% sd = 225/2;
% % μs ⟶ s
% sd = sd*1e-6; % s
% % flatness
% fl = 10;
%
% src = exp(-((t-tce)/sd).^fl);
% ------------------------------------------------------------------------------
%                                        ∂ ∫ Δt
% ------------------------------------------------------------------------------
srcd = differentiate_o6(src,dt);
srci = integrate_line(src,dt);
% ------------------------------------------------------------------------------
%                                    🎵🎵🎵🎵
% ------------------------------------------------------------------------------
[src_,f,df] = fourier_rt(src,dt);
[srcd_,f,df] = fourier_rt(srcd,dt);
[srci_,f,df] = fourier_rt(srci,dt);
% ------------------------------------------------------------------------------
%
%                                   🎨🎨🎨🎨
%
% ------------------------------------------------------------------------------
verde   = [0.4660 0.6740 0.1880];
azul    = [0.3010 0.7450 0.9330];
purpura = [0.4940 0.1840 0.5560];
naranja = [0.8500 0.3250 0.0980];
% ------------------------------------------------------------------------------
figure;
subplot(1,2,1);
hold on;
plot(t*1e6,src,'linewidth',2,'color',verde);
plot(t*1e6,src,'.','markersize',15,'color',azul);

% plot(t*1e6,srcd,'linewidth',2,'color',azul);
% plot(t*1e6,srcd,'.','markersize',15,'color',purpura);

plot(t*1e6,srci,'linewidth',2,'color',purpura);
plot(t*1e6,srci,'.','markersize',15,'color',naranja);
hold off;
axis tight;
axis square;
xlabel('Time (μs)')
ylabel('Amplitude ( )')
simple_figure()

subplot(1,2,2);
loglog(f,abs(src_)/numel(real(src_)),'linewidth',2,'color',verde);
hold on;
loglog(f,abs(srcd_)/numel(real(srcd_)),'linewidth',2,'color',azul);
loglog(f,abs(srci_)/numel(real(srci_)),'linewidth',2,'color',purpura);
hold off;
axis tight;
axis square;
xlabel('Frequency (Hz)')
ylabel('Power')
simple_figure()
% ------------------------------------------------------------------------------
