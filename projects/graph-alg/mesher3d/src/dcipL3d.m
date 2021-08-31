function V = dcipL3d(n_g2m,n_ij,n_IJ,I,J,neigh_type,graph2mesh,x,y,z,sigm,alphas)
% diego domenzain
% august 2021
% ------------------------------------------------------------------------------
alphas = alphas(x,y,z,sources_,pts_robin)
  alphas = zeros(nprobin,1);

  ipx = ;
  alpha_r = 0;
  alpha_cos=0;
  for isource=1:nsource
    isx = ;
    radi_po = sqrt( (x(isx)-x(ipx))^2 + (y(isy)-y(ipy))^2 + (z(isz)-z(ipz))^2 );
    ca_po   = abs(z(isz) - z(ipz));

    isx = ;
    radi_ne = sqrt( (x(isx)-x(ipx))^2 + (y(isy)-y(ipy))^2 + (z(isz)-z(ipz))^2 );
    ca_ne   = abs(z(isz) - z(ipz));

    alpha_r  = alpha_r + radi_ne - radi_po;
    alpha_cos= alpha_cos + (ca_ne/radi_ne^3) - (ca_po/radi_po^3);
  end
  alphas(iprobin) = alpha_r * alpha_cos;

iyxz = graph2mesh(J(il))
iyxz_= graph2mesh(J(ii))
[dx,dy,dz]= deltas(iyxz,iyxz_nei,x,y,z)
  [ix,iy,iz] = get_ixyz(iyxz,nx,ny,nz);
  [ix_,iy_,iz_] = get_ixyz(iyxz_nei,nx,ny,nz)


% diego domenzain
% August 2021 @ Aarhus University
% ------------------------------------------------------------------------------
%
% discretization of
%                   L ≈ - ∇ ⋅ sigm ∇
%
% where
%       L = sparse(I,J,V);
% ------------------------------------------------------------------------------
% NOTE:
% • 'sigm' is not taken into account right now.
% • 'dx' & 'dz' are also not taken into account (assumes dx=dz everywhere)
% • Robin boundary conditions are NOT implemented,
%   only Neumann (no-flow) are implemented.
% ------------------------------------------------------------------------------
V = zeros(n_IJ,1);
il = 1;
il_= 0;
for i_g2m=1:n_g2m
    % -- end of line
    il_ = il_ + n_ij(i_g2m)+1;
    % -- in line
    % V
    ith = 0;
    for ii=(il+1):il_
        % sum for ith entry
        ith = ith +1;
        % neighbors entries
        V(ii) = -1;
        % NOTE: the +1 and -1 should be replaced (each) with a product of:
        % 1. material properties (i.e. harmonic averages of conductivity),
        % 2. and dx_j / dz_j terms.
        %
        % to do so, use J.
        % for example, inside this loop:
        % J(il) gives the ith node,
        % J(ii) gives neighbor of ith node.
        %
        % the harmonic average would be:
        % c_ij_ = ( 2*c(J(il)) * c(J(ii)) ) / ( c(J(il)) + c(J(ii)) )
    end
    % robin nodes are summed to 'ith'
    for i_nei=1:6
        if (neigh_type(J(il),i_nei) == 0)
            ith = ith +1;
            % NOTE: the +1 should be replaced by the appropriate entry.
            % J(il) gives the ith node,
            % neigh_type(J(il),i_nei) gives neighbor of ith node.
            %
            % WARNING: be careful with the corners!!
            % 'bottom' corners are characterized by having two 0 entries in:
            %     neigh_type(J(il),1:4)
            % 'top' corners are characterized by having one 0 and one -1 in:
            %     neigh_type(J(il),1:4)
        end
    end
    % ith entry
    V(il) = ith;
    % -- begining of next line
    il = il_ + 1;
end

end
