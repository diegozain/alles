clc
clear
close all
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% set discretization param
% NOTE: if the initial guess is pretty good, do matt's trick and increase the 
% length of the signal for the scheme to converge.
nt = 4000;
dt = 2.5e-4;
t  = dt*(0:(nt-1)).';
% ------------------------------------------------------------------------------
% set harmonic observed data
fo=5; % (Hz)
% α coefficients go with cos
% β coefficients go with sin
%
% a harmonic 'uh' function looks like:
%
% uh(t) = Σj Σi αi⋅cos(2*π*t * fo*hj) + Σj Σi βi⋅sin(2*π*t * fo*hj)
%
% ------------------------------------------------------------------------------
% uo_h = 2*cos(2*pi*fo*t);

uo_h = 2*cos(2*pi*fo*t) + cos(2*pi*fo*3*t) + cos(2*pi*fo*5*t) + cos(2*pi*fo*11*t) + 2*sin(2*pi*fo*t) + sin(2*pi*fo*3*t) + sin(2*pi*fo*5*t) + sin(2*pi*fo*11*t);
% ------------------------------------------------------------------------------
% uo_s = zeros(nt,1);
uo_s = exp(-t.^2);
% ------------------------------------------------------------------------------
uo = uo_h + uo_s;
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
beta_  = 4;%0;
fo     = 5.5; % Hz
% ------------------------------------------------------------------------------
% multiples h*fo that make up the harmonics
h = [1, 3, 5, 11];
% h = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
% h = [1];
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
fprintf('--> this is how long it took me to complete\n    the forward model.\n')
% ------------------------------------------------------------------------------
[uo_,f,df] = fourier_rt(uo,dt);
[uh_,f,df] = fourier_rt(uh,dt);
% ------------------------------------------------------------------------------
figure;

subplot(1,2,1)
hold on;
plot(t,uh,'k-','linewidth',5);
plot(t,uo,'r-');
hold off;
xlabel('Time (s)')
ylabel('Voltage (V)')
title('True & initial')
simple_figure()

subplot(1,2,2)
hold on;
plot(f,abs(uh_),'k-')
plot(f,abs(uo_),'r-')
hold off;
xlabel('Frequency (Hz)')
ylabel('Power')
title('True & initial')
simple_figure()
% ------------------------------------------------------------------------------
% vis obj function
% ------------------------------------------------------------------------------
% nals_=1e2;%5e2;
% nfos_=1e2;%5e2;
% fos_   = linspace(0,6,nfos_);
% alphas_= linspace(-6,6,nals_);
% Ob_=zeros(nfos_,nals_);
% for ifos = 1:nfos_
%   for ialphas=1:nals_
%     uh     = hd_fwd(t,[alphas_(ialphas)],betas,[fos_(ifos)],h,nt_,nt__);
%     error_ = uo - uh;
%     Ob_(ifos,ialphas) = log(sum(error_.^2));
%     % Ob_(ifos,ialphas) = sum(error_.^2);
%   end
% end
% 
% figure;
% subplot(1,2,1)
% fancy_imagesc(Ob_,alphas_,fos_)
% colormap(rainbow2_cb(1))
% axis normal
% xlabel('alphas ( )')
% ylabel('fo (Hz)')
% title('Objective function')
% simple_figure()
% %{
% ------------------------------------------------------------------------------
%
% inversion
%
% ------------------------------------------------------------------------------
% -- # of iterations
niter=4;
Obs = zeros(niter,3);

ob_fos   =zeros(niter,1);
fos_niter=zeros(nb,niter);

ob_alphas=zeros(niter,1);
alphas_niter=zeros(nb*nh,niter);
ob_betas=zeros(niter,1);
betas_niter=zeros(nb*nh,niter);
% type of objective function.
% sse = sum of squared errors
% lnsse = ln( sum of squared errors )
type_obj='sse';
% debug
steps_= zeros(niter,3);
al_fo = zeros(niter+1,2);
al_fo(1,:) = [alphas(1),fos(1)];
% -- parabola step sizes hyperparameters
nparabo   = 2;

% fo
k_fos_    = -1e-8; % -1e-8
k_fos__   =  1e-6; %  1e-6

% α & β 
k_alphas_ = -1e-3;
k_alphas__=  1e-1;
k_betas_  = -1e-3;
k_betas__ =  1e-1;

% ?
k__ =-1e-6;
k___= 1e-4;
% ------------------------------------------------------------------------------
fprintf('\n\n     Ok Mr User. I will now begin the inversion,');
fprintf('\n        it will go for %i iterations,', niter);
fprintf('\n        and will use the objective function %s.\n\n',type_obj)
% ------------------------------------------------------------------------------
tic;
for iiter=1:niter
  % -- fos
  type_obj='lnsse';
  uh          = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
  [Ob,error_] = hd_obj(uo,uh,type_obj);
  
  ob_fos(iiter,1)=Ob;
  fos_niter(:,iiter)=fos;
  % fos gradient
  g_fos = hd_grad_f(error_,t,alphas,betas,fos,h,nt_,nt__);
  % fos step size
  % % uo_s = zeros
  % k_fos_=1e-8; k_fos__=1e-4; nparabo_fos=1e2;
  % % uo_s = exp(-t^2)
  k_fos_=1e-8; k_fos__=2e-5; nparabo_fos=1e2;
  step_fos = hd_step_f(uo,uh,g_fos,Ob,k_fos_,k_fos__,nparabo_fos,type_obj,t,alphas,betas,fos,h,nt_,nt__);
  % update
  d_fos = - step_fos * g_fos;
  % frequency always positive
  fos   = fos.*exp(fos.*d_fos);
end
[~,iiter] = min(ob_fos);
fos = fos_niter(:,iiter);

for iiter=1:niter
  type_obj='sse';
  uh          = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
  [Ob,error_] = hd_obj(uo,uh,type_obj);
  
  ob_alphas(iiter,1)=Ob;
  alphas_niter(:,iiter)=alphas;
  
  ob_betas(iiter,1)=Ob;
  betas_niter(:,iiter)=betas;
  
  % -- alphas
  % alphas gradient
  g_alphas = hd_grad_a(error_,t,alphas,fos,h,nt_,nt__);
  % alphas step size
  k_alphas_=1e-8; k_alphas__=1e-3; nparabo_ab=1e2;
  step_alphas = hd_step_a(uo,uh,g_alphas,Ob,k_alphas_,k_alphas__,nparabo_ab,type_obj,t,alphas,betas,fos,h,nt_,nt__);

  % -- betas
  % betas gradient
  g_betas = hd_grad_b(error_,t,betas,fos,h,nt_,nt__);
  % betas step size
  k_betas_=1e-8; k_betas__=5e-3; nparabo_ab=1e2;
  step_betas = hd_step_b(uo,uh,g_betas,Ob,k_betas_,k_betas__,nparabo_ab,type_obj,t,alphas,betas,fos,h,nt_,nt__);
  
  % -- update
  % alphas
  d_alphas = - step_alphas * g_alphas;
  alphas   = alphas + d_alphas;
  % betas
  d_betas  = - step_betas * g_betas;
  betas    = betas + d_betas;
end
[~,iiter] = min(ob_alphas);
alphas = alphas_niter(:,iiter);
[~,iiter] = min(ob_betas);
betas = betas_niter(:,iiter);
toc;
fprintf('--> this is how long it took me to complete\n    the inversion.\n')
% ------------------------------------------------------------------------------
% recovered signal
uh = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
[uh_,f,df] = fourier_rt(uh,dt);
% ------------------------------------------------------------------------------
% % debug
% hold on;
% plot(al_fo(:,1),al_fo(:,2),'k.-','markersize',20);
% hold off;

% subplot(1,2,2)
% surf(alphas_,fos_,Ob_,'edgecolor','none')
% xlabel('alphas ( )')
% ylabel('fo (Hz)')
% zlabel('Objective function')
% simple_figure()
% ------------------------------------------------------------------------------
figure;

subplot(1,2,1)
hold on;
plot(t,uh,'k-','linewidth',3);
plot(t,uo,'r-');
hold off;
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Recovered data')
simple_figure()

subplot(1,2,2)
hold on;
plot(f,abs(uh_),'k-','linewidth',3)
plot(f,abs(uo_),'r-')
hold off;
xlabel('Frequency (Hz)')
ylabel('Power')
title('Recovered data')
simple_figure()
% ------------------------------------------------------------------------------
figure;
subplot(2,1,1)
plot(fos_niter.','.-','markersize',20)
subplot(2,1,2)
plot(ob_fos,'.-','markersize',20)
%
% figure;
% subplot(2,1,1)
% plot(betas_niter.','.-','markersize',20)
% subplot(2,1,2)
% plot(ob_betas,'.-','markersize',20)
% ------------------------------------------------------------------------------
figure;
subplot(1,2,1)
hold on;
plot(t,uo_s,'r','linewidth',4);
plot(t,uo-uh,'k','linewidth',2);
plot(t,uo,'b');
hold off;
xlabel('Time (s)')
ylabel('Voltage (V)')
simple_figure()

[uh_,f,df] = fourier_rt(uo-uh,dt);
[uo_s_,f,df]= fourier_rt(uo_s,dt);

subplot(1,2,2)
hold on;
plot(f,abs(uo_s_),'r-','linewidth',4)
plot(f,abs(uh_),'k-','linewidth',2)
plot(f,abs(uo_),'b-')
hold off;
xlabel('Frequency (Hz)')
ylabel('Power')
simple_figure()
% ------------------------------------------------------------------------------
% window mean width = (1/smallest harmonic) / dt
[uo_,~] = window_mean(uo,800);
[uh_,~] = window_mean(uo-uh,800);
[uh__,~]= window_mean_(uo-uh,800);

figure;
hold on;
plot(t,uo_s,'r','linewidth',4);
plot(t,uo-uh,'b','linewidth',2);
plot(t,uh_,'k');
plot(t,uh__,'c--');
plot(t,uo_,'g--')
hold off;

figure;
hold on;
plot(uo_s,'r','linewidth',4);
plot(uo-uh,'b','linewidth',2);
plot(uh_,'k.-','markersize',20);
plot(uh__,'c.-','markersize',20);
hold off;


% [u,C]= window_mean_(ones(10,1),3);
u=[-9 5 3 -4 0 1 2 -6 0 8];
v=[1 1 1 1 1]/5;
conv(u,v,'same')
[u_,C]= window_mean_(u.',5);
u_.'

u=[-9 5 3 -4 0 1 2 -6 0 8];
v=[1 1 1 1]/4;
conv(u,v,'same')
[u_,C]= window_mean_(u.',4);
u_.'

%}
