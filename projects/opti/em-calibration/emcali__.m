clear
close all
clc
% ------------------------------------------------------------------------------
addpath('src/')
% ------------------------------------------------------------------------------
%
%
%
%                          🌐 grid-search inversion 🌐
%
%
%
%
% ------------------------------------------------------------------------------
% 📩
% load('..\..\..\..\tmp\TEM_coil_cal.mat');
load('..\..\..\..\tmp\TEM_coil_cal_NG1.mat');
% load('..\..\..\..\tmp\TEM_coil_cal_NG2.mat');
% ------------------------------------------------------------------------------
r= complex(double(R),0);
a= complex(double(Rd1),0);
b= complex(double(Rd2),0);
s= s.';
datao1 = D1.';
datao2 = D2.';
fprintf('\n\n        a=%2.2f • b=%2.2f\n',a,b);
clear R Rd1 Rd2 D1 D2
% ------------------------------------------------------------------------------
datao = datao2 ./ datao1;
ns = numel(s);
% ------------------------------------------------------------------------------
fprintf('\n\n   gonna put some noise 🎶\n\n')
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
%
%
%                                👾 inversion 👾
%
%
% ------------------------------------------------------------------------------
% 📦
param = zeros(5,1);
param = complex(param,0);
param(3) = r;
param(4) = a;
param(5) = b;
% ------------------------------------------------------------------------------
% 📟
niter = 500; % 100;
tolerr = 1e-6;
obj = tolerr + 1;
obj_ = zeros(niter,1);
steps_ = zeros(niter,1);
% ------------------------------------------------------------------------------
% 🌐
%  1e-6  < l < 1e-2
%  1e-12 < c < 1e-9
llmin = -6;
llmax = 1;
ccmin = -12;
ccmax = -6;

nll=5;
ncc=5;
obje = zeros(nll,ncc);
% ------------------------------------------------------------------------------
% 📞🎨
objeall = zeros(nll,ncc,niter);
ccall = zeros(ncc,niter);
llall = zeros(nll,niter);

parammm=zeros(2,niter);
% ------------------------------------------------------------------------------
iter_=1;
tic;
while (obj>tolerr & iter_ < niter+1)
  % 🌐🔍
  ll = logspace(llmin,llmax,nll);
  cc = logspace(ccmin,ccmax,ncc);
  obje = zeros(nll,ncc);
  for ill=1:nll
    param(1) = ll(ill);
    for icc=1:ncc
      param(2) = cc(icc);
      % 👉
      data = fwdemcali(param,s);
      % Θ & residual
      [obj,resi] = objemcali(data,datao);
      % ⬛
      obje(ill,icc) = obj;
    end
  end
  [obj,iobj] = min(obje(:));
  % 🌐 ⟵
  [ill,icc] = ind2sub([nll,ncc],iobj);
  param(1) = ll(ill);
  param(2) = cc(icc);
  llmin = log10(ll(ill) - ll(ill)*0.5/iter_);
  llmax = log10(ll(ill) + ll(ill)*0.5/iter_);
  ccmin = log10(cc(icc) - cc(icc)*0.5/iter_);
  ccmax = log10(cc(icc) + cc(icc)*0.5/iter_);
  % 🚶🚶
  obj_(iter_) = obj;
  iter_=iter_+1;
  % 🐾
  if (mod(iter_-1,100)==0)
    fprintf('    ♠ just did iteraion #%i\n',iter_-1);
  end
  % 📞🎨
  parammm(:,iter_) = param(1:2);
  objeall(:,:,iter_-1) = obje;
  ccall(:,iter_-1) = cc;
  llall(:,iter_-1) = ll;
end
toc;
% ------------------------------------------------------------------------------
fprintf('\n\n             😎😎😎😎😎😎 recovered\n l = %2.2d\n c = %2.2d\n             😎😎😎😎😎😎\n\n',param(1),param(2));
% ------------------------------------------------------------------------------
%                                  🎨🎨
% ------------------------------------------------------------------------------
verde   = [0.4660 0.6740 0.1880];
azul    = [0.3010 0.7450 0.9330];
purpura = [0.4940 0.1840 0.5560];
naranja = [0.8500 0.3250 0.0980];
amarillo= [0.9290 0.6940 0.1250];
% ------------------------------------------------------------------------------
data = fwdemcali(param,s);

figure(1);
subplot(1,2,1);
hold on;
semilogx(imag(s),real(datao),'-','linewidth',3,'color',[0.4940 0.1840 0.5560]);
semilogx(imag(s),real(data),'-','linewidth',3,'color',[0.4660 0.6740 0.1880]);
hold off;
axis tight;
axis square;
grid on;
xlabel('iω (rad)')
ylabel('R data ( - )')
simple_figure();

subplot(1,2,2);
hold on;
semilogx(imag(s),imag(datao),'-','linewidth',3,'color',[0.4940 0.1840 0.5560]);
semilogx(imag(s),imag(data),'-','linewidth',3,'color',[0.4660 0.6740 0.1880]);
hold off;
axis tight;
axis square;
grid on;
legend({'observed','recovered'})
xlabel('iω (rad)')
ylabel('iR data ( - )')
simple_figure();
% ------------------------------------------------------------------------------
figure;
semilogy(obj_,'k.-','markersize',20);
axis tight;
axis square;
xlabel('Iteration #')
ylabel('Objective function')
simple_figure()
% ------------------------------------------------------------------------------
figure;
hold on;
fancy_imagesc(obje,log10(cc),log10(ll));
for iter_=1:niter
  if (obj_(iter_)~=0)
    fancy_imagesc(objeall(:,:,iter_),log10(ccall(:,iter_)),log10(llall(:,iter_)));
  end
end
plot(log10(param(2)),log10(param(1)),'w.','markersize',40)
plot(log10(parammm(2,1)),log10(parammm(1,1)),'.','markersize',15,'color',azul)
for iter_=1:niter
  if (obj_(iter_)~=0)
    plot(log10(parammm(2,iter_)),log10(parammm(1,iter_)),'.','markersize',10,'color',[0,0,0])
    plot(log10(parammm(2,iter_)),log10(parammm(1,iter_)),'.','markersize',5,'color',naranja)
  end
end
plot(log10(param(2)),log10(param(1)),'.','markersize',20,'color',purpura)
hold off;
colormap(rainbow2_cb(1))
xlabel('lg C')
ylabel('lg L')
simple_figure()
% ------------------------------------------------------------------------------
