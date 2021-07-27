function g_betas = hd_grad_b(error_,t,betas,fos,h,nt_,nt__)

nt = numel(t);
nb = numel(fos);
nh = numel(h);

% betas gradient
error_  = error_.';
nbh     = nb*nh;
g_betas= zeros(nbh,1);

for ibh=1:nbh
  % translate one for-loop into two for-loops
  ih = mod(ibh,nh);
  if (ih==0)
    ih=nh;
  end
  ib = ((ibh - ih) / nh) + 1;

  % time-interval times
  indexes_t = (1 + (ib-1)*(nt_-nt__)):(nt_+ (ib-1)*(nt_-nt__));
  % sine block
  argu_     = 2*pi* (fos(ib)*t(indexes_t)) * h(ih);
  sin_bloc_ = sin( argu_ );
  % dot product
  g_betas(ibh) = error_(indexes_t) * sin_bloc_;
end

end