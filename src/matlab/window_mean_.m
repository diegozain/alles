function [v,C] = window_mean_(u,nw)
% diego domenzain. jul 2021 @ AU
% window-mean the columns of u,
% with a width c.
% ------------------------------------------------------------------------------
% this construction is better than conv(u,window mean,'same') because
% it takes into account averages at early and late times.
% ------------------------------------------------------------------------------
nt= numel(u);
C = zeros(nt,nt);

if mod(nw,2)==0;
  nb_ = nw/2;
  % last block starts at (-1)
  nb__ = nt - nb_;
else
  nb_ = ceil(nw/2);
  % last block starts at (-1)
  nb__ = nt - nb_ + 1;
end
% ------------------------------------------------------------------------------
%
% -- build convolution without matrix --
%
% ------------------------------------------------------------------------------
v = zeros(nt,1);

% starting block
for it=1:(nb_-1)
  v(it) = (ones(1,it+nb_-1)/(it+nb_-1)) * u(1:(it+nb_-1));
end
% middle block
for it=nb_:nb__
  v(it) = (ones(1,nw)/nw) * u((it-nb_+1):(it-nb_+nw));
end
% ending block
for it=(nb__+1):nt
  v(it) = (ones(1,(nt-it+nb_))/(nt-it+nb_)) * u((it-nb_+1):nt);
end
% ------------------------------------------------------------------------------
%
% -- build convolution matrix --
%
%    v = C*u
% ------------------------------------------------------------------------------
if nargout > 1
  % starting block
  for it=1:(nb_-1)
    C(it,1:(it+nb_-1)) = ones(1,it+nb_-1)/(it+nb_-1);
  end
  % middle block
  for it=nb_:nb__
    C(it,(it-nb_+1):(it-nb_+nw)) = ones(1,nw)/nw;
  end
  % ending block
  for it=(nb__+1):nt
    C(it,(it-nb_+1):nt) = ones(1,(nt-it+nb_))/(nt-it+nb_);
  end
end
end
