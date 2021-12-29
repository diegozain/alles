function g_fos = hd_grad_f(error_,t,alphas,betas,fos,h,nt_,nt__)
% diego domenzain. jul 2021

nt = numel(t);
nb = numel(fos);
nh = numel(h);

% fos gradient
error_  = error_.';
g_fos   = zeros(nb,1);

for ib=1:nb
  % time-interval times of size (nt_ × 1)
  indexes_t = (1 + (ib-1)*(nt_-nt__)):(nt_+ (ib-1)*(nt_-nt__));
  % harmonic interval of size (1 × nh)
  indexes_h = (1 + (ib-1)*nh):(nh + (ib-1)*nh);

  % argument of size nt_ × nh
  argu_     = 2*pi* (fos(ib)*t(indexes_t)) * h;
  % cosine block of size nt_ × nh
  cos_bloc_ = (-2*pi*t(indexes_t)*h).*sin( argu_ );
  % sine block of size nt_ × nh
  sin_bloc_ = (2*pi*t(indexes_t)*h).*cos( argu_ );
  % multiply by vec of size nh,
  % result is a vec of size nt_
  cos_bloc__= cos_bloc_ * alphas(indexes_h);
  sin_bloc__= sin_bloc_ * betas(indexes_h);

  % dot product
  g_fos(ib) = error_(indexes_t) * (cos_bloc__ + sin_bloc__);
end
end
