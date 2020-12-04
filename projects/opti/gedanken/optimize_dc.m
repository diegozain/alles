close all
clear
clc
% ..............................................................................
addpath('src/')
% ..............................................................................
% source, parameters and data
s_dc  = [1;-1];
p_dc_o= [12; 23];
Mdc   = [1 -1];
[nd_dc,nu_dc] = size(Mdc);
% ..............................................................................
[u_dc_o,d_dc_o] = fwd_dc(p_dc_o,s_dc,Mdc,zeros(nd_dc,1));
% ..............................................................................
% parameter search space
%
np = 5e+1;
p1 = linspace(3,15,np);
p2 = linspace(3,24,np);
% ..............................................................................
% objective function
%
Ndc = eye(length(d_dc_o));
% ..............................................................................
Nreg = eye(length(p_dc_o));
% ..............................................................................
EE = Edc_surface(p1,p2,s_dc,Mdc,d_dc_o,Ndc);
% ..............................................................................
% check
%
E1_ = min(EE');
E2_ = min(EE);
[~,i] = min(E1_);
[~,ii] = min(E2_);
pp = [p1(i); p2(ii)];
% ..............................................................................
fprintf('p grid search \n')
fprintf('%d \n',pp)
fprintf('\n')
fprintf('p_o \n')
fprintf('%d \n',p_dc_o)
fprintf('\n')
fprintf('data from p_o and Ldc \n')
fprintf('%d \n',d_dc_o)
fprintf('\n')
%-------------------------------------------------------------------------------
%
%                            optimization
%                               joint
%-------------------------------------------------------------------------------
p = [5; 5];
% ..............................................................................
fprintf('p seed \n')
fprintf('%d \n',p)
fprintf('\n')
% ..............................................................................
max_iter = 100;
[p,P,E,ss] = gd_dc(p,s_dc,Mdc,d_dc_o,Ndc,max_iter);
% ..............................................................................
fprintf('p optimized \n')
fprintf('%d \n',p)
fprintf('\n')
% ..............................................................................
[~,no_iter] = size(P);
% ..............................................................................
fprintf('no. iter & E \n')
fprintf('%1d %d \n',no_iter,E(end))
fprintf('\n')
%-------------------------------------------------------------------------------
figure;
imagesc(p2,p1,log10(normali(EE)))
hold on
plot(P(2,1:end),P(1,1:end),'.-','Markersize',15)
plot(p(2),p(1),'.y','Markersize',40)
plot(P(2,1),P(1,1),'.b','Markersize',60)
plot(p_dc_o(2),p_dc_o(1),'pr','Markersize',30,'MarkerFaceColor','red')
hold off
colormap(flipud(bone))
% caxis([0, 1e-2])
frame = gca;
frame.YTickLabel = [];
frame.XTickLabel = [];
frame.XTick = [];
frame.YTick = [];
colorbar off
ylabel('Parameter 2')
xlabel('Parameter 1')
title('Objective function # 1')
simple_figure()
% ------------------------------------------------------------------------------
% fig = gcf;fig.PaperPositionMode = 'auto'; print(gcf,'optim-dc','-dpng','-r600')
% ------------------------------------------------------------------------------
% ..............................................................................
%{
% ..............................................................................
figure;
hold on
plot(P(2,1:end),P(1,1:end),'.-','Markersize',10)
plot(p(2),p(1),'.y','Markersize',40)
plot(P(2,1),P(1,1),'.b','Markersize',60)
plot(p_dc_o(2),p_dc_o(1),'pr','Markersize',30,'MarkerFaceColor','red')
hold off
% xlim([1 15])
% ylim([1 15])
set(gca,'Ydir','reverse')
ylabel('p_1')
xlabel('p_2')
title('p search path')
% ..............................................................................
figure;
hold on
plot(1,E(1),'.b','Markersize',40)
for i=2:10:no_iter-1
    plot(i,E(i),'.','Markersize',20)
end
plot(no_iter,E(end),'.y','Markersize',40)
hold off
ylabel('E')
xlabel('iteration #')
title('error history')
% ..............................................................................
figure;
hold on
plot(1,ss(1),'.b','Markersize',40)
for i=2:10:no_iter-1
    plot(i,ss(i),'.','Markersize',20)
end
plot(no_iter,ss(end),'.y','Markersize',40)
hold off
ylabel('norm( g )')
xlabel('iteration #')
title('gradient norm history')
% ..............................................................................
%}
% ..............................................................................
