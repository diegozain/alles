function g_alphas = hd_grad_a(error_,t,alphas,fos,h,nt_,nt__)

nt = numel(t);
nb = numel(fos);
nh = numel(h);

% alphas gradient
error_  = error_.';
nbh     = nb*nh;
g_alphas= zeros(nbh,1);

for ibh=1:nbh
 % translate one for-loop into two for-loops
 ih = mod(ibh,nh);
 if (ih==0)
   ih=nh;
 end
 ib = ((ibh - ih) / nh) + 1;

 % time-interval times
 indexes_t = (1 + (ib-1)*(nt_-nt__)):(nt_+ (ib-1)*(nt_-nt__));
 % cosine block
 argu_     = 2*pi* (fos(ib)*t(indexes_t)) * h(ih);
 cos_bloc_ = cos( argu_ );
 % dot product
 g_alphas(ibh) = error_(indexes_t) * cos_bloc_;
end
end