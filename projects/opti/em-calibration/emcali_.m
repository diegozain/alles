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
load('C:\Users\au699065\OneDrive - Aarhus Universitet\Documents\diego\code\tmp\TEM_coil_cal.mat');
% ------------------------------------------------------------------------------
r= complex(double(R),0);
a= complex(double(Rd1),0);
b= complex(double(Rd2),0);
s= s.';
datao1 = D1.';
datao2 = D2.';

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
niter = 200; % 100;
% 🚚 true values
%          l = 1e-3;
%          c = 1e-10;
%  1e-6  < l < 1e-2
%  1e-12 < c < 1e-9
l= 1e-1; % 1e-2
c= 1e-8; % 1e-9
kparam_ = 1e-7;  % 1e-7;
kparam__= 1e-1;  % 1e-1;
nparabo = 3;     % 3
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
fprintf('\n\n             😎😎😎😎😎😎\n l = %2.2d\n c = %2.2d\n             😎😎😎😎😎😎\n\n',param(1),param(2));
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

% subplot(1,3,3)
% hold on;
% plot(real(data),imag(data),'-','linewidth',3,'color',[0.9290 0.6940 0.1250]);
% hold off;
% axis tight;
% % axis square;
% simple_figure();
% ------------------------------------------------------------------------------
obj_ = zeros(niter,1);
steps_ = zeros(niter,1);
% ------------------------------------------------------------------------------
dparam_=zeros(5,1);
iter_=1;
while (iter_ < niter+1)
  % 👉
  data = fwdemcali(param,s);
  % Θ & residual
  [obj,resi] = objemcali(data,datao);
  % ∇
  grad_ = grademcali(param,s,resi);
  % step-size
  step_ = stepemcali(obj,param,datao,grad_,s,kparam_,kparam__,nparabo);
  % p ⟵ p + Δp
  % dparam = -step_*grad_ + dparam_;
  dparam = -step_*grad_;
  param = param.*exp(param.*dparam);
  % 🚶🚶
  % dparam_=1e-1*dparam;
  obj_(iter_) = obj;
  steps_(iter_) = step_;
  iter_=iter_+1;
end
% ------------------------------------------------------------------------------
fprintf('\n\n             😎😎😎😎😎😎\n l = %2.2d\n c = %2.2d\n             😎😎😎😎😎😎\n\n',param(1),param(2));
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

% subplot(1,3,3)
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
ll = logspace(-6,3,nll);
cc = logspace(-12,-8,ncc);
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
hold off;
axis normal;
colormap(rainbow2_cb(1))
xlabel('lg C')
ylabel('lg L')
simple_figure()
% ------------------------------------------------------------------------------
