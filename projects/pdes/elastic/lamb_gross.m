function [vx_rx,vz_rx] = lamb_anal_v3(rx,sz,Fo,vp,vs,rho,f,u)

% rx = dx*(1:50).';
% vp=vp(1);vs=vs(1);rho=rho(1);
% [sz,f,df] = fourier_rt_(fz_,dt);
% u = linspace(0,1/(vel_min*0.4),100).';
% [vx_rx,vz_rx] = lamb_anal_v3(rx,f,sz,Fo,vp,vs,rho,u);
% figure;plot(real(vz_rx(:,2)))

mu = rho*vs^2;

% sz = fftshift(sz);

nrx= numel(rx);
nf = numel(f);
nu = numel(u);

w = 2*pi*f;
wrx = rx * w;
wrx = wrx.';

vx_rx = zeros(nf,nrx);
vz_rx = zeros(nf,nrx);

du=u(2)-u(1);

% ------------------------------------------------------------------------------
tuk= tukey(2*nu,2);
tuk= tuk((nu+1):(2*nu));
u  = u.*tuk;

for iu=1:nu
  u_ = u(iu);
  
  uwrx = u_*wrx;
  
  a = sqrt( (1/vp^2) - u_^2 );
  b = sqrt( (1/vs^2) - u_^2 );
  
  u_vx = (u_*(u_^2 - b^2 + 2*a*b) / ((u_^2 - b^2)^2 + 4*a*b*u_^2)) * sin(uwrx);
  u_vz = (1i*u_*a / (vs^2 * (u_^2 - b^2)^2 + 4*a*b*u_^2)) * cos(uwrx);
  
  vx_rx = vx_rx + u_vx*du;
  vz_rx = vz_rx + u_vz*du;
end

wpimu = (1i*w/(pi*mu)).';
wpimu = repmat(wpimu,1,nrx);

vx_rx = Fo * wpimu .* vx_rx;
vz_rx = - Fo * wpimu .* vz_rx;

vx_rx = ifft(vx_rx,[],1);
vz_rx = ifft(vz_rx,[],1);


vx_rx=real(vx_rx);
vz_rx=real(vz_rx);
end