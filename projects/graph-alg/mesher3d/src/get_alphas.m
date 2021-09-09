function alphas = get_alphas(x,y,z,srcs_xyz,robin_xyz)
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
% compute α to get the right boundary conditions at subsurface nodes
% when solving:
%                      -∇ ⋅ σ ∇ ϕ = s
%
% this function 'get_alphas' assumes 'srcs_xyz' is a set of srcs_xyz done at the
% same moment (i.e. an 'ab' pair, but also supports multi-source schemes).
% ------------------------------------------------------------------------------
% srcs_xyz  : (nsources) × (xyz) × (±)
% robin_xyz : (nprobin) × (xyz)
%
% both of these contain indexes in the mesh, NOT the graph.
% ------------------------------------------------------------------------------
nprobin= size(robin_xyz,1);
alphas = zeros(nprobin,1);
for iprobin=1:nprobin
  ipx = robin_xyz(iprobin,1);
  ipy = robin_xyz(iprobin,2);
  ipz = robin_xyz(iprobin,3);

  alpha_r = 0;
  alpha_cos=0;
  for isource=1:nsource
    % ⚪ positive source
    isx = srcs_xyz(isource,1,1);
    isy = srcs_xyz(isource,2,1);
    isz = srcs_xyz(isource,3,1);
    radi_po = sqrt( (x(isx)-x(ipx))^2 + (y(isy)-y(ipy))^2 + (z(isz)-z(ipz))^2 );
    ca_po   = abs(z(isz) - z(ipz));
    % ⚫ negative source
    isx = srcs_xyz(isource,1,2);
    isy = srcs_xyz(isource,2,2);
    isz = srcs_xyz(isource,3,2);
    radi_ne = sqrt( (x(isx)-x(ipx))^2 + (y(isy)-y(ipy))^2 + (z(isz)-z(ipz))^2 );
    ca_ne   = abs(z(isz) - z(ipz));
    % α = Σ ∓ r(s±) ⋅ Σ ∓ ca(s±) / r(s±)^3
    alpha_r  = alpha_r + radi_ne - radi_po;
    alpha_cos= alpha_cos + (ca_ne/radi_ne^3) - (ca_po/radi_po^3);
  end
  alphas(iprobin) = alpha_r * alpha_cos;
end
end
