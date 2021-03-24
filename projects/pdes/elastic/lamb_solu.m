function [vx_rx,vz_rx] = lamb_solu(rx,sz,Fo,vp,vs,rho,f,u)
% diego domenzain
% @ CSM, spring 2021
% ------------------------------------------------------------------------------
% taken from lamb.f:
% https://git.scc.kit.edu/Seitosh/Seitosh/-/blob/master/src/synt/misc/lamb.f
%
% 
%
% ------------------------------------------------------------------------------
% rx:       receiver distances from source (in m)
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
% rx = dx*(10:100).';
% vp=vp(1);vs=vs(1);rho=rho(1);
% [sz,f,df] = fourier_rt_(fz_,dt);
% % u = linspace(0,1.2e-3,10000).'; 
% u = linspace(0,1/(vel_min),100).';
% [vx_rx,vz_rx] = lamb_solu(rx,sz,Fo,vp,vs,rho,f,u);
% 
% figure; fancy_imagesc(vx_rx); axis normal
% figure; fancy_imagesc(vz_rx); axis normal
% ------------------------------------------------------------------------------
nrx = numel(rx);
nf = numel(f);
nu = numel(u);

vx_rx = zeros(nf,nrx);
vz_rx = zeros(nf,nrx);
% ------------------------------------------------------------------------------
mu = rho*vs^2;
du = u(2) - u(1);

% - rayleigh thing
raylim = 1e-30;
% ------------------------------------------------------------------------------
tuk  = tukey(2*nu,0.9);
tuk  = tuk((nu+1):(2*nu));
utap = u.*tuk;

tuk = tukey(nf,1);
tuk = tuk.';
ftuk= f .* fftshift(tuk);
% ftuk=f;
% ------------------------------------------------------------------------------
for iu=1:nu
  u_ = utap(iu);
  utap_=utap(iu);
  
  a = sqrt( (1/vp^2) - u_^2 );
  b = sqrt( (1/vs^2) - u_^2 );
  % ----------------------------------------------------------------------------
  % taper
  
  % ----------------------------------------------------------------------------
  % rayleigh ??
  rayleigh = (2*u_^2 - 1/vs^2)^2 + (4*u_^2);
  % ----------------------------------------------------------------------------
  % calculate the numerators for the greens function
  vx_num = (u_^2) * ((1/vs^2) - 2*u_^2 - 2*a*b);
  vz_num = u_ * (1/vs^2) * a * 1i;
  % ----------------------------------------------------------------------------
  % frequency loop
  omegau = 2*pi*f*u_;
  jarg=rx*omegau;
  jarg=jarg.';
  
  J0=besselj(0,jarg);
  J1=besselj(1,jarg);
  
  Ivx = J1 * ((vx_num/rayleigh)*utap_*du);
  Ivz = J0 * ((vz_num/rayleigh)*utap_*du);
  
  vx_rx=vx_rx + Ivx;
  vz_rx=vz_rx + Ivz;
end
% vx_rx_ = ifft(vx_rx,[],1);
% vx_rx_ = real(vx_rx_);
% figure; fancy_imagesc(vx_rx_); axis normal
% figure;plot(vx_rx_(:,10))
% ------------------------------------------------------------------------------
% sources
src = ftuk.*f*(-Fo/mu) .* sz.';
src = src.';
src = repmat(src,1,nrx);

vx_rx=vx_rx.*src;
vz_rx=vz_rx.*src;
% ------------------------------------------------------------------------------
% ifft
vx_rx = ifft(vx_rx,[],1);
vz_rx = ifft(vz_rx,[],1);

vx_rx= - real(vx_rx);
vz_rx= - real(vz_rx);
end % function lamb_solu
