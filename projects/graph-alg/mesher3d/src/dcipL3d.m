function V = dcipL3d(n_g2m,n_ij,n_IJ,I,J,neigh_type,graph2mesh,robin_graph,alphas,n_ar,sig,x,y,z)
% diego domenzain
% august 2021
% ------------------------------------------------------------------------------
%
% discretization of
%                   L ≈ - ∇ ⋅ σ ∇
%
% where
%       L = sparse(I,J,V);
% ------------------------------------------------------------------------------
%
% the ith row of matrix L has the following form,
%
%       L(i,:) = [ δik·σi + Σj (δij · σ⋆ij)        (-δij · σ⋆ij) ]
%                         ith entry                  jth entries
%
% 'j' runs over all inner neighbors of node 'i'.
% 'k' runs over all boundary neighbors of node 'i'.
%
%           σ⋆ij = 2·σi·σj / (σi + σj)              harmonic average
%            δij = Δij / Δ⟂ij                   inner & air-ground nodes
%            δik = Δki·αi                subsurface boundary nodes (0 otherwise)
%            δik = Δk1·c1 + Δk2·c2              corner nodes (0 otherwise)
% ------------------------------------------------------------------------------
% n_g2m      • total # of nodes in the graph.
% n_ij       • holds the info of how many entries in I (and J) belong to each
%              node i. It is an array of size n_g2m by 1.
% n_IJ       • total number of non-zero entries of L (also length of I and J).
% I          • I denotes the row entries
% J          • J the column entries
% neigh_type •
% graph2mesh •
% mesh2graph •
% x          •
% y          •
% z          •
% sigm       •
% alphas     •
% ------------------------------------------------------------------------------
% NOTE:
% • 'sigm' is not taken into account right now.
% • 'dx' & 'dz' are also not taken into account (assumes dx=dz everywhere)
% • Robin boundary conditions are NOT implemented,
%   only Neumann (no-flow) are implemented.
% ------------------------------------------------------------------------------
% dij  = deltas(J(il),J(ij),x,y,z);
%      iyxz = graph2mesh(J(il))
%      iyxz_= graph2mesh(J(ij))
%      [ix,iy,iz]    = get_ixyz(iyxz,nx,ny,nz)
%      [ix_,iy_,iz_] = get_ixyz(iyxz_,nx,ny,nz)
% ------------------------------------------------------------------------------
%    L(i_g2m,:) = [ δik·σi + Σj (δij · σ⋆ij)        (-δij · σ⋆ij) ]
%                         ith entry                  jth entries
% ------------------------------------------------------------------------------
V = zeros(n_IJ,1);
% ------------------------------------------------------------------------------
% begining of row in L
il = 1;
% end of row in L (initialized)
il_= 0;

% beginning of robins
iar = 1;
% end of robins (initialized)
iar_= 0;
iprobin_ = 1;
% ------------------------------------------------------------------------------
% ⚫ access a node 'i_g2m' in the graph,
%    which is indexed by 'il' in I, J, & V.
for i_g2m=1:n_g2m
    % -- end of row in L
    il_ = il_ + n_ij(i_g2m)+1;
    % -- in row of L
    % 'ith' will be the entry of L(i_g2m,i_g2m).
    %  here the value is reset to zero.
    ith = 0;
    % ◼️ loop thru inner neighbors of 'i_g2m',
    %    whose position in I, J, & V is indexed by 'ij',
    %    → asign values to,
    %
    %                    L(i_g2m , ij) = -δij · σ⋆ij
    %
    %    → and do the inner neighbor sum for the entry L(i_g2m,i_g2m),
    %                                                            Σj (δij · σ⋆ij)
    %    inside this loop:
    %      • J(il) gives the ith node in the graph,
    %      • J(ij) gives inner neighbor of ith node in the graph.
    % --
    for ij=(il+1):il_
        % δij : neighbor faces and edges of FV scheme.
        % dij = deltas(J(il),J(ij),x,y,z);
        dij = 1;
        % σ⋆ij : harmonic average
        sig_ij_ = ( 2*sig(J(il)) * sig(J(ij)) ) / ( sig(J(il)) + sig(J(ij)) );
        % sum for the entry L(i_g2m , i_g2m)
        ith = ith + dij*sig_ij_;
        % entry of L(i_g2m , ij)
        V(ij) = -dij*sig_ij_;
    end
    % ◻️ loop thru robin neighbors of 'i_g2m',
    %    whos position in I, J, & V does not exist!
    %    to find them, you need to access:
    %                       neigh_type(J(il),i_nei)
    %    if the number you find is 0,
    %    then position i_nei is robin for 'i_g2m'.
    %
    %    → gives the second part of the entry L(i_g2m,i_g2m),
    %                                                    δik·σi
    %    inside this loop:
    %      • J(il) gives the ith node in the graph.
    %
    % ⚠️ robin nodes are summed to 'ith' for now.
    % --
    % is the current node a robin node?
    if (J(il) == robin_graph( iar ,1))
      % -- end of robin node
      iar_ = iar_ + n_ar( iprobin_ );
      for iiar=iar:iar_
        % -- get αi
        alphai = alphas( iiar );
        % % -- get δik
        % % type of neighbor for this robin node,
        % %  1 2 3 4 5  6
        % %  → ↑ ← ↓ ⦿ ⊛
        % i_nei = robin_graph( iiar ,2);
        % dik = deltas_robin(J(il),i_nei,x,y,z);
        dik = 1;

        ith = ith + dik*alphai;
      end
      % beginning of next robin node
      iar = iar_ + 1;
      iprobin_ = iprobin_ + 1;
    end
    % ⚫ now that we have all the information about the ith entry,
    %    i.e. L(i_g2m,i_g2m), we can plug it in.
    %
    %             L(i_g2m,i_g2m) = δik·σi + Σj (δij · σ⋆ij)
    % --
    V(il) = ith;
    % -- begining of next row in L
    il = il_ + 1;
end

end
