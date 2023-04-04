clc
clear
close all
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% set discretization param
nt = 8000;
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
% more complicated harmonic
uo_h = 2*cos(2*pi*fo*t) + cos(2*pi*fo*3*t) + cos(2*pi*fo*5*t) + cos(2*pi*fo*11*t) + 2*sin(2*pi*fo*t) + sin(2*pi*fo*3*t) + sin(2*pi*fo*5*t) + sin(2*pi*fo*11*t);
% ------------------------------------------------------------------------------
% actual signal
% uo_s = zeros(nt,1);
% uo_s = t.^2;
uo_s = exp(-t.^2) + 5;
% uo_s = cos(2*pi*fo*3.5*t);
% ------------------------------------------------------------------------------
% observed signal
uo   = uo_h + uo_s;
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
% multiples h*fo that make up the harmonics
h = [1, 3, 5, 11];
nh= numel(h);
% ------------------------------------------------------------------------------
%                      initial α, β, & fo
% ------------------------------------------------------------------------------
% alpha_ = 4;
% beta_  = 4;
% fo     = 5.5; % Hz
%
% alphas = ones(nb*nh,1) * alpha_;
% betas  = ones(nb*nh,1) * beta_;
% fos    = ones(nb,1) * fo;
foi = fo;
% foi = fo*0.99;
% h=[1 3];nh=2;
% ------------------------------------------------------------------------------
% initial guess for α & β
alpha_ = 1*std(uo);
beta_  = 1*std(uo);
% ----------------------------------------------------------------------------
alphas = 1*ones(nb*nh,1) * alpha_ ./ (1:nb*nh).';
betas  = 1*ones(nb*nh,1) * beta_  ./ (1:nb*nh).';
fos    = ones(nb,1) * foi;
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
rgb=qualitcolor(4);

figure;

subplot(1,2,1)
hold on;
plot(t,uh,'linewidth',2,'color',rgb(1,:));
plot(t,uo,'color',rgb(2,:));
hold off;
xlabel('Time (s)')
ylabel('Voltage (V)')
legend({'initial harmonic','observed data'})
axis tight;
axis square;
simple_figure()

subplot(1,2,2)
loglog(f,abs(uh_),'linewidth',2,'color',rgb(1,:))
hold on;
loglog(f,abs(uo_),'color',rgb(2,:))
hold off;
xlabel('Frequency (Hz)')
ylabel('Power')
axis tight;
axis square;
xticks([1,10,100,1000])
yticks([1,10,100,1000])
grid on;
simple_figure()
% ------------------------------------------------------------------------------
%
% 👈 inversion
%
% ------------------------------------------------------------------------------
% -- # of iterations
niter=5;
% ------------------------------------------------------------------------------
fprintf('\n\n     Ok Mr User. I will now begin the inversion,');
fprintf('\n        it will go for %i iterations.\n\n', niter);
% ------------------------------------------------------------------------------
tic;
[alphas,betas,fos] = hd_inversion(uo,niter,t,alphas,betas,fos,h,nt_,nt__);
toc;
fprintf('--> this is how long it took me to complete\n    the inversion.\n')
% ------------------------------------------------------------------------------
% -- recovered signal
% harmonic
uh = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
[uh_,f,df] = fourier_rt(uh,dt);
% signal (without window mean, i.e. without makeup)
us = uo-uh;
% ------------------------------------------------------------------------------
% 
%                                   🖌️🎨
% 
% ------------------------------------------------------------------------------
figure;

subplot(1,2,1)
hold on;
plot(t,uh,'linewidth',2,'color',rgb(1,:));
plot(t,uo,'color',rgb(2,:));
hold off;
xlabel('Time (s)')
ylabel('Voltage (V)')
legend({'recovered harmonic','observed data'})
axis tight;
axis square;
simple_figure()

subplot(1,2,2)
loglog(f,abs(uh_),'linewidth',2,'color',rgb(1,:))
hold on;
loglog(f,abs(uo_),'color',rgb(2,:))
hold off;
xlabel('Frequency (Hz)')
ylabel('Power')
axis tight;
axis square;
xticks([1,10,100,1000])
yticks([1,10,100,1000])
grid on;
simple_figure()
% ------------------------------------------------------------------------------
figure;
subplot(1,2,1)
hold on;
plot(t,uo_s,'linewidth',4,'color',rgb(1,:));
plot(t,us,'linewidth',2,'color',rgb(2,:));
plot(t,uo,'color',rgb(3,:));
hold off;
xlabel('Time (s)')
ylabel('Voltage (V)')
legend({'true signal','recovered data','observed data'})
axis tight;
axis square;
simple_figure()

[uo_s_,f,df]= fourier_rt(uo_s,dt);
[uh_,f,df] = fourier_rt(us,dt);

subplot(1,2,2)
loglog(f,abs(uo_s_),'linewidth',4,'color',rgb(1,:))
hold on;
loglog(f,abs(uh_),'linewidth',2,'color',rgb(2,:))
plot(f,abs(uo_),'color',rgb(3,:))
hold off;
xlabel('Frequency (Hz)')
ylabel('Power')
axis tight;
axis square;
xticks([1,10,100,1000])
yticks([1,10,100,1000])
grid on;
simple_figure()
% ------------------------------------------------------------------------------
% window mean width = (1/smallest harmonic) / dt
us_ = window_mean_(us,ceil(1/min(fos)/dt));
% ------------------------------------------------------------------------------
figure;

subplot(1,2,1)
hold on;
plot(t,uo_s,'linewidth',4,'color',rgb(1,:));
plot(t,us,'linewidth',2,'color',rgb(2,:));
plot(t,us_,'--','linewidth',2,'color',rgb(3,:));
plot(t,uo,'linewidth',1,'color',rgb(4,:));
hold off;
xlabel('Time (s)')
ylabel('Voltage (V)')
axis tight;
axis square;
simple_figure()

[us__,f,df]= fourier_rt(us_,dt);

subplot(1,2,2)
loglog(f,abs(uo_s_),'linewidth',4,'color',rgb(1,:))
hold on;
loglog(f,abs(uh_),'linewidth',2,'color',rgb(2,:))
loglog(f,abs(us__),'--','linewidth',2,'color',rgb(3,:))
loglog(f,abs(uo_),'color',rgb(4,:))
hold off;
xlabel('Frequency (Hz)')
ylabel('Power')
legend({'true signal','recovered data','conv. rec. data','observed data'})
axis tight;
axis square;
xticks([1,10,100,1000])
yticks([1,10,100,1000])
grid on;
simple_figure()
% ------------------------------------------------------------------------------
%}
