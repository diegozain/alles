close all
clear
clc
% ------------------------------------------------------------------------------
fprintf('\n     ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■')
fprintf('\n            some text here')
fprintf('\n     ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■\n')
% ------------------------------------------------------------------------------
t = linspace(0,2*pi,10000);
dt= t(2)-t(1);
y = sin(t);
y_= differentiate_o6(y,dt);
% ------------------------------------------------------------------------------
figure('visible', 'off');
hold on;
plot(t,y,'linewidth',2)
plot(t,y_,'linewidth',2)
hold off;
axis tight;
axis square;
xlabel('Time')
ylabel('Amplitude')
simple_figure()
% ------------------------------------------------------------------------------
print(gcf,'../../pics/somepic','-dpng','-r350')
% ------------------------------------------------------------------------------
fprintf('\n     ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■')
fprintf('\n          some text here too')
fprintf('\n     ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■ ■\n')
% ------------------------------------------------------------------------------
% remember to close all figures before exiting so matlab 👶 doesnt cry 😢
close all;
quit;
% ------------------------------------------------------------------------------
