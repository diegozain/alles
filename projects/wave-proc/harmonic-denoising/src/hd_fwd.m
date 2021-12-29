function uh = hd_fwd(t,alphas,betas,fos,h,nt_,nt__)
% diego domenzain

nt = numel(t);
nb = numel(fos);
nh = numel(h);

% harmonic signal
uh = zeros(nt,1);
for ib=1:nb
  % --- one block ---
  % each block is of size nt_ x nh
  % cos_bloc = cos( 2*pi*fo*t*h )

  % time-interval times
  indexes_t = (1 + (ib-1)*(nt_-nt__)):(nt_+ (ib-1)*(nt_-nt__));
  % harmonic interval
  indexes_h = (1 + (ib-1)*nh):(nh + (ib-1)*nh);

  % this is actually a matrix of size nt_ x nh
  argu_ = 2*pi* (fos(ib)*t(indexes_t)) * h;

  cos_bloc = cos( argu_ );
  sin_bloc = sin( argu_ );

  uh_=zeros(nt,1);
  uh_(indexes_t) = cos_bloc * alphas(indexes_h) + sin_bloc * betas(indexes_h);

  % record the harmonic signal WITH overlapping blocks
  uh = uh + uh_;
end
end
