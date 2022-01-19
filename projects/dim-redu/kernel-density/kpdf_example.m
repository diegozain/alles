clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% example from wikipedia
% https://en.wikipedia.org/wiki/Kernel_density_estimation
datao = [-2.1; -1.3; -0.4; 1.9; 5.1; 6.2];
% ------------------------------------------------------------------------------
% [kpdf,x] = kernel_density(datao);
[kpdf,x] = kernel_density(datao,2.25);
% ------------------------------------------------------------------------------
figure;
hold on;
plot(x,kpdf,'linewidth',3);
plot(datao,zeros(numel(datao)),'.','markersize',30)
hold off;
axis tight;
xlabel('Data space');
ylabel('Density function');
simple_figure();
% ------------------------------------------------------------------------------
