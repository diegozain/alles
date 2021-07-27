function [u,C] = window_mean_(u,nw)
% diego domenzain. jul 2021 @ AU
% window-mean the columns of u,
% with a width c.
% ---
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

u=C*u;
end