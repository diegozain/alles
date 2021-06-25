function V = div_c_grad(n_g2m,n_ij,n_IJ,neigh_type,I,J,c)
% diego domenzain
% @ Colorado School of Mines
% spring 2021
% ------------------------------------------------------------------------------
% 
% discretization of 
%                   L ≈ - ∇ ⋅ c ∇
% 
% where 
%       L = sparse(I,J,V);
% ------------------------------------------------------------------------------
% NOTE:
% • 'c' is not taken into account right now.
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
    for i_nei=1:4
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