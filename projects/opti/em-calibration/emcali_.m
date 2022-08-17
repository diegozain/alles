clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src/')
% ------------------------------------------------------------------------------
% d = ( l*c*s.^2 + (r*c + l/a)*s + 1 + r/a ) ./ ( l*c*s.^2 + (r*c + l/b)*s + 1 + r/b )
%
%                           Θ = || do - d ||_2^2
%                           ∇ = (∂l , ∂c)
%
% ∂l(d) = (b*s*(b - a)) ./ (a*(b*c*s*(l*s + r) + b + l*s + r)^2)
% ∂c(d) = (b*s*(a - b).*(l*s + r)^2) ./ (a*(b*c*s.*(l*s + r) + b + l*s + r)^2)
% ------------------------------------------------------------------------------
% 📩
% load('..\..\..\..\tmp\TEM_coil_cal.mat');
% load('..\..\..\..\tmp\TEM_coil_cal_NG1.mat');
load('..\..\..\..\tmp\TEM_coil_cal_NG2.mat');
% ------------------------------------------------------------------------------
r= complex(double(R),0);
a= complex(double(Rd1),0);
b= complex(double(Rd2),0);
s= s.';
datao1 = D1.';
datao2 = D2.';
fprintf('\n\n        a=%2.2f • b=%2.2f\n\n',a,b);
clear R Rd1 Rd2 D1 D2
% ------------------------------------------------------------------------------
datao = datao2 ./ datao1;
ns = numel(s);
% ------------------------------------------------------------------------------
fprintf('\n\n   gonna put some noise 🎶\n')
rng(1);
noimagre = std(real(datao));
noimagim = std(imag(datao));
noimagre = 0.1*noimagre;
noimagim = 0.1*noimagim;
datao = datao + (noimagre*rand(numel(s),1) + noimagim*1i*rand(numel(s),1));
% ------------------------------------------------------------------------------
% fprintf('\n\n   gonna halve the data 💕\n')
% datao = datao(1:2:ns);
% s = s(1:2:ns);
% ------------------------------------------------------------------------------
% 📟
niter = 500; % 100;
% 🚚 true values
%          l = 1e-3;
%          c = 1e-10;
%  1e-6  < l < 1e-2
%  1e-12 < c < 1e-9
l= 10^(-4.4); % 1e-5; % 8e-5; % 1e-1; % 1e-2
c= 10^(-8.6); % 1e-7; % 1.6e-9; % 1e-8; % 1e-9
kparam_ = 1/(ns*1e4); % 1e-9;  % 1e-7;
kparam__= 1e2/ns; % 1e-3;  % 1e-1;
nparabo = 3;     % 3
% ------------------------------------------------------------------------------
% % the idea here is that given a value for 'l',
% % then 'c' has to be a certain way,
% % so maybe doing what is beneath makes the initial guess better.
% % not really though 😞
% ------------------------------------------------------------------------------
% for ii=1:10
%   c = (r/a - 1 + l*s.*(1/a - datao/b)) ./ ((datao-1)*l.*s.^2 + r*s.*(datao - 1));
%   c = abs(real(c(ns)));
%   l = (r/a - 1 + r*c*s.*(1 - datao)) ./ ((datao-1)*c.*s.^2 + s.*(datao/b - 1/a));
%   l = abs(real(l(ns)));
% end
% ------------------------------------------------------------------------------
% 📦
param = zeros(5,1);
param = complex(param,0);
param(1) = l;
param(2) = c;
param(3) = r;
param(4) = a;
param(5) = b;

paramm=param;
% ------------------------------------------------------------------------------
fprintf('\n\n             😎😎😎😎😎😎 initial\n l = %2.2d\n c = %2.2d\n             😎😎😎😎😎😎\n\n',param(1),param(2));
% ------------------------------------------------------------------------------
data = fwdemcali(param,s);

figure(1);
subplot(1,2,1);
semilogx(imag(s),real(data),'-','linewidth',3,'color',[0.9290 0.6940 0.1250]);
axis tight;
axis square;
xlabel('iω (rad)')
ylabel('R data ( - )')
simple_figure();

subplot(1,2,2);
semilogx(imag(s),imag(data),'-','linewidth',3,'color',[0.9290 0.6940 0.1250]);
axis tight;
axis square;
xlabel('iω (rad)')
ylabel('iR data ( - )')
simple_figure();
% ------------------------------------------------------------------------------
% figure;
% hold on;
% plot(real(data),imag(data),'-','linewidth',3,'color',[0.9290 0.6940 0.1250]);
% hold off;
% axis tight;
% axis square;
% xlabel('R data ( - )')
% ylabel('iR data ( - )')
% simple_figure();
% ------------------------------------------------------------------------------
tolerr = 1e-5;
obj = tolerr + 1;
obj_ = zeros(niter,1);
steps_ = zeros(niter,1);

parammm=zeros(2,niter);
% ------------------------------------------------------------------------------
iter_=1;
tic;
while (obj>tolerr & iter_ < niter+1)
  % 👉
  data = fwdemcali(param,s);
  % Θ & residual
  [obj,resi] = objemcali(data,datao);
  % ∇
  grad_ = grademcali(param,s,resi);
  % 📐
  step_ = stepemcali(obj,param,datao,grad_,s,kparam_,kparam__,nparabo);
  % p ⟵ p + Δp
  dparam = -step_*grad_;
  % ☝️📅
  param = param.*exp(param.*dparam);

  % 🚶🚶
  obj_(iter_) = obj;
  steps_(iter_) = step_;
  iter_=iter_+1;
  % 📞
  parammm(:,iter_-1) = param(1:2);
end
toc;
% ------------------------------------------------------------------------------
fprintf('\n\n             😎😎😎😎😎😎 recovered\n l = %2.2d\n c = %2.2d\n             😎😎😎😎😎😎\n\n',param(1),param(2));
% ------------------------------------------------------------------------------
data = fwdemcali(param,s);

figure(1);
subplot(1,2,1);
hold on;
semilogx(imag(s),real(data),'-','linewidth',3,'color',[0.4660 0.6740 0.1880]);
semilogx(imag(s),real(datao),'-','linewidth',3,'color',[0.4940 0.1840 0.5560]);
hold off;
axis tight;
axis square;
grid on;
xlabel('iω (rad)')
ylabel('R data ( - )')
simple_figure();

subplot(1,2,2);
hold on;
semilogx(imag(s),imag(data),'-','linewidth',3,'color',[0.4660 0.6740 0.1880]);
semilogx(imag(s),imag(datao),'-','linewidth',3,'color',[0.4940 0.1840 0.5560]);
hold off;
axis tight;
axis square;
grid on;
legend({'initial','recovered','observed'})
xlabel('iω (rad)')
ylabel('iR data ( - )')
simple_figure();
% ------------------------------------------------------------------------------
% figure(2);
% hold on;
% plot(real(data),imag(data),'-','linewidth',3,'color',[0.4660 0.6740 0.1880]);
% plot(real(datao),imag(datao),'--','linewidth',3,'color',[0.4940 0.1840 0.5560]);
% hold off;
% axis tight;
% axis square;
% legend({'initial','recovered','observed'})
% xlabel('R data ( - )')
% ylabel('iR data ( - )')
% simple_figure();
% ------------------------------------------------------------------------------
figure;
semilogy(obj_,'k.-','markersize',20);
axis tight;
axis square;
xlabel('Iteration #')
ylabel('Objective function')
simple_figure()
% ------------------------------------------------------------------------------
% l= 1e-3;
% c= 1e-10;
%  1e-6  < l < 1e-2
%  1e-12 < c < 1e-9
nll=1e2;
ncc=1e2;
% ll = logspace(-4,-2,nll);
% cc = logspace(-11,-9,ncc);
ll = logspace(-6,1,nll);
cc = logspace(-12,-6,ncc);
obje = zeros(nll,ncc);
param_=param;
for ill=1:nll
  param_(1) = ll(ill);
  for icc=1:ncc
    param_(2) = cc(icc);
    % 👉
    data = fwdemcali(param_,s);
    % Θ & residual
    [obj,resi] = objemcali(data,datao);
    % ⬛
    obje(ill,icc) = obj;
  end
end
figure;
hold on;
fancy_imagesc(obje,log10(cc),log10(ll));
plot(log10(param(2)),log10(param(1)),'w.','markersize',40)
plot(log10(paramm(2)),log10(paramm(1)),'c.','markersize',40)
for iter_=1:niter
  plot(log10(parammm(2,iter_)),log10(parammm(1,iter_)),'c.','markersize',5)
end
hold off;
colormap(rainbow2_cb(1))
xlabel('lg C')
ylabel('lg L')
simple_figure()
% ------------------------------------------------------------------------------
% 🐛
% lichao gave me 3 examples with the same 'l' & 'c',
% but different 'a' & 'b' and he didn't give me a heads up. fucker.
% so this plot was made to figure that out.
figure(4);
subplot(1,2,1);
hold on;
semilogx(imag(s),real(datao),'-','linewidth',3);
hold off;
axis tight;
axis square;
grid on;
xlabel('iω (rad)')
ylabel('iR data ( - )')
simple_figure();

subplot(1,2,2);
hold on;
semilogx(imag(s),imag(datao),'-','linewidth',3);
hold off;
axis tight;
axis square;
grid on;
xlabel('iω (rad)')
ylabel('iR data ( - )')
simple_figure();
% ------------------------------------------------------------------------------
