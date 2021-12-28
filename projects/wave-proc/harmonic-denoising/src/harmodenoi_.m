function us = harmodenoi_(uo,dt,fo,h)
% diego domenzain
% ago 2021
% ------------------------------------------------------------------------------
% given an observed signal 'uo' with harmonic noise, find the signal 'us'.
%
% uo :: signal
% dt :: time sampling
% fo :: frequency
% h  :: numbers for which fo gets its harmonics
% ------------------------------------------------------------------------------
% α coefficients go with cos
% β coefficients go with sin
%
% a harmonic 'uh' function looks like:
%
% uh(t) = Σj Σi αi⋅cos(2*π*t * fo*hj) + Σj Σi βi⋅sin(2*π*t * fo*hj)
% ------------------------------------------------------------------------------
nt= numel(uo);
nh= numel(h);

t= (0:1:(nt-1))*dt;
t=t.';
% ------------------------------------------------------------------------------
% nb  : # of time blocks
% nt_ : number of samples in the over-lapping time
% nt__: overlapping number of time samples   (nt__ = ceil(1/fo/dt))
nb = 1;
[nb,nt_,nt__] = hd_nbnt_(nt,fo,dt,nb);
% ------------------------------------------------------------------------------
% initial guess for α & β
alpha_ = 1*std(uo);
beta_  = 1*std(uo);
% ------------------------------------------------------------------------------
alphas = ones(nb*nh,1) * alpha_ ./ (1:nb*nh).';
betas  = ones(nb*nh,1) * beta_  ./ (1:nb*nh).';
fos    = ones(nb,1) * fo;
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
%
% inversion
%
% ------------------------------------------------------------------------------
% -- # of iterations
% between 4 and 10 seems to work ok for all examples so far
niter=6;
[alphas,betas,fos] = hd_inversion(uo,niter,t,alphas,betas,fos,h,nt_,nt__);
% ------------------------------------------------------------------------------
% -- recovered signal
% harmonic
uh = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
% signal (without window mean, i.e. without makeup)
us = uo-uh;
% ------------------------------------------------------------------------------
% window mean width = (1/smallest harmonic) / dt
% nw = fix((1/(min(fos)*h(1)))/dt);
nw = fix((1/(fo*h(1)))/dt);
us = window_mean_(us,nw);
end
