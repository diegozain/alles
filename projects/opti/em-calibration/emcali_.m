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
% ------------------------------------------------------------------------------
% 📟
niter = 1000; % 1000;
% true values
% l= 1e-3;
% c= 1e-10;
%  1e-6  < l < 1e-2
%  1e-12 < c < 1e-9
l= 1e-2;
c= 1e-9;
kparam_ = 1e-10; % 1e-10;
kparam__= 1e-2; % 1e-2;
nparabo = 5; % 3
% ------------------------------------------------------------------------------
% 📦
param = zeros(5,1);
param = complex(param,0);
param(1) = l;
param(2) = c;
param(3) = r;
param(4) = a;
param(5) = b;
% ------------------------------------------------------------------------------
obj_ = zeros(niter,1);
% ------------------------------------------------------------------------------
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
  param = param.*exp(- step_*param.*grad_);
  % 🚶🚶
  obj_(iter_) = obj;
  iter_=iter_+1;
end
% ------------------------------------------------------------------------------
fprintf('\n\n             😎😎😎😎😎😎\n l = %2.2d\n c = %2.2d\n             😎😎😎😎😎😎\n\n',param(1),param(2));
% ------------------------------------------------------------------------------
figure;
semilogy(obj_,'k.-','markersize',20);
axis tight;
axis square;
xlabel('Iteration #')
ylabel('Objective function')
simple_figure()
% ------------------------------------------------------------------------------
