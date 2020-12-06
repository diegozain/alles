function b = beamformer_(fo,r,d_,f,v,theta)
% diego domenzain
% fall 2018 @ BSU
% ..............................................................................
% b = zeros(ntheta,nv);
[nf,nr] = size(d_);
nv = numel(v);
ntheta = numel(theta);
k = [cos(theta).' sin(theta).'];

% kv*r.' is of size (theta by r)
% kv*r.' * d_(ifo,:).' is of size (theta by 1)
% 
b = zeros(ntheta,nv);
for iv=1:nv
  kv = k/v(iv);
  b(:,iv) = exp(1i*2*pi*fo*kv*r.') * (d_(binning(f,fo),:)');
end
b = b/nr;
end