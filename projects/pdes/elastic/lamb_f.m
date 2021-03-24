function [vx_rx,vz_rx] = lamb_f(rx,sz,F,vp,vs,rho,f,u)
% diego domenzain
% @ CSM, spring 2021
% ------------------------------------------------------------------------------
% taken from lamb.f:
% https://git.scc.kit.edu/Seitosh/Seitosh/-/blob/master/src/synt/misc/lamb.f
%
% 
%
% ------------------------------------------------------------------------------
% r:        receiver distances from source (in m)
% sz:       fourier transform of source on vz
% F:        strength of the force pointing into the halfspace (in N)
% vp:       p-velocity (in m/s)
% vs:       s-velocity (in m/s)
% rho:      density (in kg/(m**3))
% f:        frequency
% u:        horizontal slowness
% ------------------------------------------------------------------------------
% a:        vertical p-slowness corresponding to u
% b:        vertical s-slowness corresponding to u
% mu:       rigidity
% vx:       horizontal displacement pointing away from source (in m)
% vz:       vertical displacement pointing into halfspace (in m)
% ------------------------------------------------------------------------------
% using my code, 
%
% r = dx*(1:50).';
% [sz,f,df] = fourier_rt_(fz_,dt);
% u = linspace(0,1.2e-3,10000).'; 
% u = linspace(0,1/(vel_min*1.4),100).';
% [vx_anal,vz_anal] = lamb_anal(r,sz,Fo,vp(1),vs(1),rho(1),f,u);
% 
% figure; fancy_imagesc(vx_anal,r,t); axis normal
% figure; fancy_imagesc(vz_anal,r,t); axis normal
% ------------------------------------------------------------------------------
nrx = numel(rx);
nf = numel(f);
nu = numel(u);

vx_rx = zeros(nf,nrx);
vz_rx = zeros(nf,nrx);
% ------------------------------------------------------------------------------
mu = rho*vs^2;

du = u(2) - u(1);
umin = u(1);
umax = u(nu);

fmin = min(abs(f));
fmax = max(abs(f));
% ------------------------------------------------------------------------------
% constants
% - cosine taper
uwil = 1e-11;
uwir = 9e-4;
% - cosine taper
fwil = 0.85;
fwir = 0.95;

% % - cosine taper
% uwil = 1e-11;
% uwir = 1e-11;
% % - cosine taper
% fwil = 1e-11;
% fwir = 1e-11;

% - rayleigh thing
raylim = 1e-30;
% ------------------------------------------------------------------------------
for iu=1:nu
  u_ = u(iu);
  
  a = sqrt( (1/vp^2) - u_^2 );
  b = sqrt( (1/vs^2) - u_^2 );
  
  ddu = du;
  if (iu==1) | (iu==nu)
    ddu=du*0.5;
  end
  % ----------------------------------------------------------------------------
  % taper
  utap = costap(umin,uwil,uwir,umax,u_);
  ddu=ddu*utap;
  % ----------------------------------------------------------------------------
  % rayleigh ??
  rayleigh = ((2*u_^2 - 1/vs^2)^2) + (4*u_^2);
  reray = real(rayleigh);
  imray = imag(rayleigh);
  absray= abs(rayleigh);
  % beware of the rayleigh-pole
  if (absray==0)
    rayleigh=raylim;
  elseif (absray < raylim)
    rayleigh= (reray/abs(reray)) + 1i*(imray/abs(imray))*raylim;
  end
  % ----------------------------------------------------------------------------
  % calculate the numerators for the greens function
  vx_num=(u_^2) * ((1/vs^2) - 2*u_^2 - 2*a*b);
  vz_num=u_ * (1/vs^2) * a * 1i;
  % ----------------------------------------------------------------------------
  % frequency loop
  for if_=1:nf
    omegau = 2*pi*f(if_)*u_;
    % --------------------------------------------------------------------------
    % receiver loop 
    for ir=1:nrx
      jarg=omegau*rx(ir);
      
      J0=besselj(0,jarg);
      J1=besselj(1,jarg);

      Iq=(J1*vx_num/rayleigh)*ddu;
      Iw=(J0*vz_num/rayleigh)*ddu;
      
      vx_rx(if_,ir)=vx_rx(if_,ir) + Iq;
      vz_rx(if_,ir)=vz_rx(if_,ir) + Iw;
    end
  end
end

% ------------------------------------------------------------------------------
% sources
for if_=1:nf
  ftap= costap(fmin,fwil,fwir,fmax,abs(f(if_)));
  ftap= sign(f(if_)) * ftap;
  src = (ftap*(-F)*f(if_)/mu) * sz(if_);
  for ir=1:nrx
    vx_rx(if_,ir)=vx_rx(if_,ir)*src;
    vz_rx(if_,ir)=vz_rx(if_,ir)*src;
  end
end
% ------------------------------------------------------------------------------
% ifft
vx_rx = ifft(vx_rx,[],1);
vz_rx = ifft(vz_rx,[],1);

vx_rx=real(vx_rx);
vz_rx=real(vz_rx);
end % function lamb_anal

% ------------------------------------------------------------------------------

% function costap     (it's a simple cosine-taper)
function costap_ = costap(min_,wil,wir,max_,val)
  if (val < min_)
    costap_=0;
  elseif (val <= wil)
    if (wil == min_)
      costap_=0;
    else
      costap_=0.5 - 0.5 * cos((val-min_)*pi/(wil-min_));
    end
  elseif (val > max_)
    costap_=0;
  elseif (val >= wir)
    if (max_ == wir)
      costap_=0;
    else
      costap_= cos((val-wir) * pi / (max_-wir));
    end
  else
    costap_=1;
  end
costap_=costap_*val;
end