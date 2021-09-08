function alphas = get_alphas(x,y,z,sources,pts_robin)
% diego domenzain
% ?.2021
% ------------------------------------------------------------------------------
% α = Σ ∓ r(s±) ⋅ Σ ∓ ca(s±) / r(s±)^3
%
%
%          .--------------------.
%         /|              🌳   /|
%        / | s.               / |
%       /  |       🏃        /  |
%      .--------------------.   |
%      |   |                |   |
%      |   . -------------- |---.
%   z  |  /                 |   /
%      | /           p.     |  / y
%      |/   🐙              |/
%      .--------------------.
%                  x
%
% ------------------------------------------------------------------------------
% sources   : (nsources) × (xyz) × (±)
% pts_robin : (nprobin) × (xyz)
% ------------------------------------------------------------------------------
nprobin= size(pts_robin,1);
alphas = zeros(nprobin,1);
for iprobin=1:nprobin
  ipx = pts_robin(iprobin,1);
  ipy = pts_robin(iprobin,2);
  ipz = pts_robin(iprobin,3);

  alpha_r = 0;
  alpha_cos=0;
  for isource=1:nsource
    % ⚪ positive source
    isx = sources(isource,1,1);
    isy = sources(isource,2,1);
    isz = sources(isource,3,1);
    radi_po = sqrt( (x(isx)-x(ipx))^2 + (y(isy)-y(ipy))^2 + (z(isz)-z(ipz))^2 );
    ca_po   = abs(z(isz) - z(ipz));
    % ⚫ negative source
    isx = sources(isource,1,2);
    isy = sources(isource,2,2);
    isz = sources(isource,3,2);
    radi_ne = sqrt( (x(isx)-x(ipx))^2 + (y(isy)-y(ipy))^2 + (z(isz)-z(ipz))^2 );
    ca_ne   = abs(z(isz) - z(ipz));
    % α = Σ ∓ r(s±) ⋅ Σ ∓ ca(s±) / r(s±)^3
    alpha_r  = alpha_r + radi_ne - radi_po;
    alpha_cos= alpha_cos + (ca_ne/radi_ne^3) - (ca_po/radi_po^3);
  end
  alphas(iprobin) = alpha_r * alpha_cos;
end
end
