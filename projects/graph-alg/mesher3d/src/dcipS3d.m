function V = dcipS3d(n_g2m,n_ij,n_IJ,I,J,neigh_mesh,graph2mesh,robin_graph,alphas,n_ar,sig,phi,x,y,z)
% diego domenzain
% august 2021
% ------------------------------------------------------------------------------
% σ 👈 ϕ,s
%
% discretization of
%
%                   S ≈ - ((∇σ L) ϕ)⊤
%
%
% ------------------------------------------------------------------------------
%
% the ith row of matrix L has the following form,
%
%   S(i,:) = [ -δik·ϕi + Σj (ϕj - ϕi)·δij·∂i(σ⋆ij)      (ϕi - ϕj)·δji·∂i(σ⋆ij) ]
%                         ith entry                          jth entries
%
% 'j' runs over all inner neighbors of node 'i'.
% 'k' runs over all boundary neighbors of node 'i'.
%
%       ∂i(σ⋆ij) = 2·σj^2 / (σi + σj)^2              harmonic average
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
% neigh_mesh • row indexes are graph nodes 🍇. Row entries are neighbors of
%              that node, in the mesh 🎲.
% graph2mesh • indexes are graph nodes 🍇, entries are mesh nodes 🎲
% robin_graph• robin nodes in the mesh 🍇. of size nprobin × 2.
%              in the second column you find which side is robin.
%              for example, a corner node will be repeated 3 times in robin_mesh
%              and the second column will perhaps read 3,4,5, meaning
%              left ←, down ↓, front ⦿ neighbors are robin.
% alphas     • alpha coefficients along robin boundary nodes
% n_ar       • counts how many entries are repeated in robin_mesh.
% sig        • σ in the 🍇
% x          • x discretization 🎲
% y          • y discretization 🎲
% z          • z discretization 🎲
% ------------------------------------------------------------------------------
% S(i,:) = [ -δik·ϕi + Σj (ϕj - ϕi)·δij·∂i(σ⋆ij)      (ϕi - ϕj)·δji·∂i(σ⋆ij) ]
%                       ith entry                          jth entries
%
% ∂i(σ⋆ij) = 2·σj^2 / (σi + σj)^2
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
    % 'ith' will be the entry of S(i_g2m,i_g2m).
    %  here the value is reset to zero.
    ith = 0;
    % ◼ loop thru inner neighbors of 'i_g2m',
    %    whose position in I, J, & V is indexed by 'ij',
    %    → asign values to,
    %
    %                    S(i_g2m , ij) = (ϕi - ϕj)·δji·∂i(σ⋆ij)
    %
    %    → and do the inner neighbor sum for the entry S(i_g2m,i_g2m),
    %                                                  Σj (ϕj - ϕi)·δij·∂i(σ⋆ij)
    %    inside this loop:
    %      • J(il) gives the ith node in the graph,
    %      • J(ij) gives inner neighbor of ith node in the graph.
    % --
    for ij=(il+1):il_
        % δij : neighbor faces and edges of FV scheme.
        iyxz = graph2mesh(J(il));
        iyxz_= graph2mesh(J(ij));
        dij  = deltas3d(J(il),neigh_mesh,iyxz,iyxz_,x,y,z);

        % ∂i(σ⋆ij) = 2·σj^2 / (σi + σj)^2
        sig_ij_ = 2*(sig(J(ij)))^2 / ( sig(J(il)) + sig(J(ij)) )^2;

        % sum for the entry S(i_g2m , i_g2m)
        ith = ith + dij*sig_ij_*(phi(J(ij)) - phi(J(il)));
        % entry of S(i_g2m , ij)
        V(ij) = dij*sig_ij_*(phi(J(il)) - phi(J(ij)));
    end
    % ◻ loop thru robin neighbors of 'i_g2m',
    %    whose position in I, J, & V does not exist!
    %    to find them, you need to access:
    %                       neigh_type(J(il),i_nei)
    %    if the number you find is 0,
    %    then position i_nei is robin for 'i_g2m'.
    %
    %    → gives the second part of the entry S(i_g2m,i_g2m),
    %                                                    -δik·ϕi
    %    inside this loop:
    %      • J(il) gives the ith node in the graph.
    %
    % ⚠ robin nodes are summed to 'ith' for now.
    % --
    % is the current node a robin node?
    if (J(il) == robin_graph( iar ,1))
      % -- end of robin node
      iar_ = iar_ + n_ar( iprobin_ );
      for iiar=iar:iar_
        % -- get αi
        alphai = alphas( iiar );
        % -- get δik
        % type of neighbor for this robin node,
        %  1 2 3 4 5  6
        %  → ↑ ← ↓ ⦿ ⊛
        i_nei = robin_graph( iiar ,2);
        iyxz = graph2mesh(J(il));
        dik = deltas_robin3d(iyxz,i_nei,x,y,z);

        ith = ith - dik*alphai*phi(J(il));
      end
      % beginning of next robin node
      iar = iar_ + 1;
      iprobin_ = iprobin_ + 1;
    end
    % ⚫ now that we have all the information about the ith entry,
    %    i.e. S(i_g2m,i_g2m), we can plug it in.
    %
    %             S(i_g2m,i_g2m) = -δik·ϕi + Σj (ϕj - ϕi)·δij·∂i(σ⋆ij)
    % --
    V(il) = ith;
    % -- begining of next row in L
    il = il_ + 1;
end

end
