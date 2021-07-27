clc
clear
close all
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
load('dcip_data.mat')

d_raw = dcip_data.data;         % V
d = dcip_data.dataNoDrift(:,3); % V
t = dcip_data.t.';              % s
dt= t(2)-t(1);

d_raw = d_raw(4001:8000);
t     = t(1:4000);

nt= numel(t);
% ------------------------------------------------------------------------------
uo = d_raw;

fo=50;
fo=5;
uo = cos(2*pi*fo*t);% + cos(2*pi*fo*3*t) + cos(2*pi*fo*5*t) + cos(2*pi*fo*11*t);
% ------------------------------------------------------------------------------
figure;
plot(t,uo);
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Raw data')
simple_figure()

[uo_,f,df] = fourier_rt(uo,dt);

figure;
plot(f,abs(uo_))
xlabel('Frequency (Hz)')
ylabel('Power')
title('Raw data')
simple_figure()
% ------------------------------------------------------------------------------
% 
% the forward model is,
%
% uh = cos_blocs * α + sin_blocs * β
%
% the cos_blocs matrix (nt x nb*nh) looks a little like this:
%
%             nh
%             |
%          _____________________________
%         |        |                    |
% nt_ ->  |    *   |_________     0     |
%         |________| <-nt__  |          |
%         |        |    *    |          | * α
%         |        |_________|          |
%         |   0                  etc    |
%         |_____________________________|
% 
% 
% each block is of size nt_ x nh.
% they all overlap on nt__ samples.
% this big matrix is of size nt x nb*nh.
% there are nb = (nt-nt__)/(nt_-nt__) blocks.
% α is of size nb*nh x 1.
% ------------------------------------------------------------------------------
% --
% nt_  :: size(1) of blocks
% nt_ = 876 for nt=16000, dt=2.5e-4, nt__=80.
% nt_ = 864 for nt=4000,  dt=2.5e-4, nt__=80.
nt_ = nt; 
% --
% nt__ :: size(1) of overlap 
% take dt * nt__ = 1/fo
%
%     --> nt__ = 1/fo/dt
nt__= 0;
% --
% nb   :: # of blocks
% each block is of size nt_ x nh
% nt = (# of blocks)*(size(1) of blocks) - collisions 
% nt = nb*nt_ - nt__*(nb - 1)
%
% -->    nb = (nt-nt__)/(nt_-nt__)
%
% nb should always be an integer!
% 
% assuming nt__ is fixed (because nt__ depends on the frequency to be found),
% we can change nt_ (always an integer) until we get an integer for nb.
% we have,
%
% nt_ = ((nt-nt__) / nb) + nt__
%
% so nb must divide nt-nt__.
%
% we can factor nt-nt__ and then look for an nb which satisfies our estimate 
% length for dt_ = nt_*dt,
%
% dt_ ≈ (((nt-nt__)/(some combination of factor(nt-nt__) )) + nt__) * dt
%
% nb = the best combination of factor(nt-nt__) that is close to dt_.
%
% whatever that combination may be, it may change the initial value for nt_!
%
% nt_ = ((nt-nt__) / nb) + nt__.
%
nb = (nt-nt__)/(nt_-nt__);
% ------------------------------------------------------------------------------
% % NOTE: nb should always be an integer. idk a clever way to ensure that :(
% nt_=920:950;
% for it=1:numel(nt_)
%   nb = (nt-nt__)/(nt_(it)-nt__);
%   fprintf('\n nt_ = %i',nt_(it));
%   fprintf('\n nb  = %2.8f\n',nb);
% end
% %{
% ------------------------------------------------------------------------------
% initial guess for α, β, and fo
alpha_ = 4;
beta_  = 0;
fo     = 5.5; % Hz
% ------------------------------------------------------------------------------
% multiples h*fo that make up the harmonics
% h = [1, 3, 5, 11];
h = [1];
nh= numel(h);
% ------------------------------------------------------------------------------
alphas = ones(nb*nh,1) * alpha_;
betas  = ones(nb*nh,1) * beta_;
fos    = ones(nb,1) * fo;
% ------------------------------------------------------------------------------
% alphas = (1./h.^2).';
% alphas = normali(alphas) * alpha_;
% alphas = repmat(alphas,nb,1);
% 
% betas = alphas;
% ------------------------------------------------------------------------------
fprintf('\n total time    = %2.2f (s)',nt*dt)
fprintf('\n interval time = %2.2f (s)',nt_*dt)
fprintf('\n overlap time  = %2.2f (s)',nt__*dt)
fprintf('\n\n # of blocks   = %i\n\n',nb)
% ------------------------------------------------------------------------------
% fwd model
% ------------------------------------------------------------------------------
tic;
% harmonic signal
uh = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
toc;
% ------------------------------------------------------------------------------
uh = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
u_ = uo - uh;

figure;
hold on;
plot(t,uh,'k-');
plot(t,uo,'r-');
hold off;
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Initial model')
simple_figure()

[u_pow_,f,df] = fourier_rt(uh,dt);

figure;
hold on;
plot(f,abs(u_pow_),'k-')
plot(f,abs(uo_),'r-')
hold off;
xlabel('Frequency (Hz)')
ylabel('Power')
title('Initial model')
simple_figure()
% ------------------------------------------------------------------------------
% vis obj function
% debug mostly
nals_=1e2;%5e2;
nfos_=1e2;%5e2;
fos_   = linspace(0,6,nfos_);
alphas_= linspace(-6,6,nals_);
Ob_=zeros(nfos_,nals_);
for ifos = 1:nfos_
  for ialphas=1:nals_
    uh     = hd_fwd(t,[alphas_(ialphas)],betas,[fos_(ifos)],h,nt_,nt__);
    error_ = uo - uh;
    Ob_(ifos,ialphas) = log(sum(error_.^2));
  end
end

figure;
fancy_imagesc(Ob_,alphas_,fos_)
axis normal
xlabel('alphas')
ylabel('fo')
simple_figure()
% %{
% ------------------------------------------------------------------------------
%
% inversion
%
% ------------------------------------------------------------------------------
% -- # of iterations
niter=30;
Obs = zeros(niter,3);
% debug
steps_= zeros(niter,3);
al_fo = zeros(niter+1,2);
al_fo(1,:) = [alphas(1),fos(1)];
% -- parabola step sizes
nparabo   = 2;
k_alphas_ = -1e-6;
k_alphas__=  1e-4;
k_fos_    = -1e-8;
k_fos__   =  1e-6;
k_betas_  = -1e-6;
k_betas__ =  1e-4;
% ------------------------------------------------------------------------------
for iiter=1:niter
  
  % -- alphas
  uh     = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
  error_ = uo - uh;
  Ob     = sum(error_.^2);
  % Ob   = -sum(log(error_.^2));
  % error_ = -1./error_;
  Obs(iiter,1)= Ob;
  % alphas gradient
  g_alphas    = hd_grad_a(error_,t,alphas,fos,h,nt_,nt__);
  % alphas step size
  step_alphas = hd_step_a(uo,g_alphas,Ob,k_alphas_,k_alphas__,nparabo,t,alphas,betas,fos,h,nt_,nt__);
  % update
  d_alphas = - 1e-1*step_alphas * g_alphas;
  alphas   = alphas + d_alphas;
  
  % % -- betas
  % uh         = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
  % error_     = uo - uh;
  % Ob         = norm(error_);
  % Obs(iiter,2)= Ob;
  % % betas gradient
  % g_betas    = hd_grad_b(error_,t,betas,fos,h,nt_,nt__);
  % % betas step size
  % step_betas = hd_step_b(error_,g_betas,k_betas,t,alphas,betas,fos,h,nt_,nt__);
  % % update
  % d_betas    = - step_betas * g_betas;
  % betas      = betas + d_betas;
  
  % -- fos
  uh     = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
  error_ = uo - uh;
  Ob     = sum(error_.^2);
  % Ob    = -sum(log(error_.^2));
  % error_= -1./error_;
  Obs(iiter,3)= Ob;
  % fos gradient
  g_fos = hd_grad_f(error_,t,alphas,betas,fos,h,nt_,nt__);
  % fos step size
  step_fos = hd_step_f(uo,g_fos,Ob,k_fos_,k_fos__,nparabo,t,alphas,betas,fos,h,nt_,nt__);
  % update
  d_fos = - 1e-1*step_fos * g_fos;
  % frequency always positive
  fos   = fos.*exp(fos.*d_fos);

  % -- debug
  steps_(iiter,1)=step_fos;
  steps_(iiter,2)=step_alphas;
  % steps_(iiter,3)=step_betas;
  
  al_fo(iiter+1,:) = [alphas(1),fos(1)];
end
% ------------------------------------------------------------------------------
% debug
hold on;
plot(al_fo(:,1),al_fo(:,2),'k.-','markersize',20);
hold off;

figure;
figure;surf(alphas_,fos_,Ob_,'edgecolor','none')
% ------------------------------------------------------------------------------
uh = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
u_ = uo - uh;

figure;
hold on;
plot(t,uh,'k-','linewidth',5);
plot(t,uo,'r-');
hold off;
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Recovered data')
simple_figure()

[u_pow_,f,df] = fourier_rt(uh,dt);

figure;
hold on;
plot(f,abs(u_pow_),'k-','linewidth',5)
plot(f,abs(uo_),'r-')
hold off;
xlabel('Frequency (Hz)')
ylabel('Power')
title('Recovered data')
simple_figure()

figure;plot(Obs)
% ------------------------------------------------------------------------------
%}
