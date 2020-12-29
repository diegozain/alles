% ------------------------------------------------------------------------------
% 
% optimize a portfolio assuming markowitz model
%
%  
% ------------------------------------------------------------------------------
clear
close all
clc
% ..............................................................................
addpath('src/')
% ------------------------------------------------------------------------------
% some stock return prices.
% dimensions are time (end of month) by stock (IBM WMT SEHI)
d = [-0.091   0.019  -0.118;
0.318   0.069   0.067;
-0.107   -0.118  -0.063;
-0.037   0.010   0.533;
0.197   0.025   0.183;
-0.028   0.000   0.494;
0.015   -0.056  -0.059;
-0.073   0.145   0.305;
-0.049   -0.140  -0.045;
-0.082   0.032  -0.362;
0.178   0.038  -0.079;
0.071   0.073   0.029];
% ------------------------------------------------------------------------------
% expected return
rho_ = 50;
% how much should the weights sum to?
miu_ = 1;
% ------------------------------------------------------------------------------
[nt,nstock] = size(d);
t=1:nt;

figure;
subplot(121)
plot(t,d)
axis tight
xlabel('Time (month)')
ylabel('Stock value')
simple_figure()
% ------------------------------------------------------------------------------
S = cov(d);
r = mean(d,1);
r = r.';
ns= size(d,2);
e_= ones(ns,1);
% ------------------------------------------------------------------------------
[w,betas] = markowitz__(d,rho_,miu_);

E1=w.'*S*w;
E2=r.'*w-rho_;
E3=e_.'*w-miu_;
fprintf('  using classic Lagrange multipliers, the objective functions give:\n   risk (should be small)      %2.2d \n   return (should be zero)     %2.2d \n   sum to one (should be zero) %2.2d\n\n',E1,E2,E3)
% ------------------------------------------------------------------------------
subplot(122)
hold on
for is=1:ns
 plot(is,w(is),'.','markersize',40)
end
hold off
axis tight
xlabel('Asset #')
ylabel('Weight value')
simple_figure()
% ------------------------------------------------------------------------------
%{
% total number of iterations
niter= 500;
[w,W,E,steps_] = markowitz_(d,rho_,miu_,niter);

E1=w.'*S*w;
E2=r.'*w-rho_;
E3=e_.'*w-miu_;
fprintf('  using classic Lagrange multipliers, the objective functions give:\n   risk (should be small)      %2.2d \n   return (should be zero)     %2.2d \n   sum to one (should be zero) %2.2d\n\n',E1,E2,E3)
% ------------------------------------------------------------------------------
figure;
plot(1:numel(E),log10(E),'.-','markersize',10)
axis tight
xlabel('Iteration #')
ylabel('log of Objective function')
simple_figure()
% 
% figure;
% plot(1:numel(steps_),steps_,'.-','markersize',10)
% axis tight
% xlabel('Iteration #')
% ylabel('Step sizes')
% simple_figure()

figure;
plot(1:size(W,2),W.','.-','markersize',10)
axis tight
xlabel('Iteration #')
ylabel('Value of weights')
simple_figure()
%}
% ------------------------------------------------------------------------------
