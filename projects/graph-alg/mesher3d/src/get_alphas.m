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
% srcs_xyz  : (nsources) × (xyz) × (±) . indexes in the mesh cube 🎲
% robin_xyz : (nprobin) × (xyz)        . indexes in the mesh cube 🎲
% ------------------------------------------------------------------------------
% neighbor types. these are in the 4rth column of robin_xyz.
%      2  6
%      | /
% 3 -- i -- 1
%    / |
%   5  4
% the numbers 1, 2, 3, 4, 5 and 6,
% represent neighbors right, up, left, down, front and back.
% ------------------------------------------------------------------------------
nprobin= size(robin_xyz,1);
nsource= size(srcs_xyz,1);

alphas = zeros(nprobin,1);
for iprobin=1:nprobin
  irobx = robin_xyz(iprobin,1);
  iroby = robin_xyz(iprobin,2);
  irobz = robin_xyz(iprobin,3);
  % type of robin node.
  irob_ = robin_xyz(iprobin,4);

  alpha_r  = 0;
  alpha_cos= 0;
  for isource=1:nsource
    % ---
    % ⚪ positive source
    isx = srcs_xyz(isource,1,1);
    isy = srcs_xyz(isource,2,1);
    isz = srcs_xyz(isource,3,1);
    radi_po = sqrt( (x(isx)-x(irobx))^2 + (y(isy)-y(iroby))^2 + (z(isz)-z(irobz))^2 );

    % adjacent side depends on the kind of neighbor.
    if (irob_==1 || irob_==3 || irob_==5 || irob_==6)
      ca_po = sqrt( (x(isx)-x(irobx))^2 + (y(isy)-y(iroby))^2 );
    end
    if (irob_==2 || irob_==4)
      ca_po = abs(z(isz) - z(irobz));
    end
    % ---
    % ⚫ negative source
    isx = srcs_xyz(isource,1,2);
    isy = srcs_xyz(isource,2,2);
    isz = srcs_xyz(isource,3,2);
    radi_ne = sqrt( (x(isx)-x(irobx))^2 + (y(isy)-y(iroby))^2 + (z(isz)-z(irobz))^2 );

    % adjacent side depends on the kind of neighbor.
    if (irob_==1 || irob_==3 || irob_==5 || irob_==6)
      ca_ne = sqrt( (x(isx)-x(irobx))^2 + (y(isy)-y(iroby))^2 );
    end
    if (irob_==2 || irob_==4)
      ca_ne = abs(z(isz) - z(irobz));
    end
    % ---
    % α = Σ ∓ r(s±) ⋅ Σ ∓ ca(s±) / r(s±)^3
    alpha_r  = alpha_r + radi_ne - radi_po;
    alpha_cos= alpha_cos + (ca_ne/radi_ne^3) - (ca_po/radi_po^3);
  end
  alphas(iprobin) = alpha_r * alpha_cos;
end
end
