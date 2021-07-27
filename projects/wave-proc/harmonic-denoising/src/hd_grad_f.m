function g_fos = hd_grad_f(error_,t,alphas,betas,fos,h,nt_,nt__)

nt = numel(t);
nb = numel(fos);
nh = numel(h);

% fos gradient
error_  = error_.';
g_fos   = zeros(nb,1);

for ib=1:nb
  % time-interval times
  indexes_t = (1 + (ib-1)*(nt_-nt__)):(nt_+ (ib-1)*(nt_-nt__));
  % harmonic interval
  indexes_h = (1 + (ib-1)*nh):(nh + (ib-1)*nh);
  
  % argument
  argu_     = 2*pi* (fos(ib)*t(indexes_t)) * h;
  % cosine block
  cos_bloc_ = (-2*pi*t(indexes_t)*h).*sin( argu_ );
  cos_bloc__= cos_bloc_ * alphas(indexes_h);
  % sine block
  sin_bloc_ = (2*pi*t(indexes_t)*h).*cos( argu_ );
  sin_bloc__= sin_bloc_ * betas(indexes_h);
  
  % dot product
  g_fos(ib) = error_(indexes_t) * (cos_bloc__ + sin_bloc__);
end
end