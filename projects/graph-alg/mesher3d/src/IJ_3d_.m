function [I,J] = IJ_3d_(n_g2m,n_ij,n_IJ,neigh_graph)
% diego domenzain
% August 2021 @ Aarhus University
%
% ------------------------------------------------------------------------------
% let L be a PDE operator acting on the graph of a grid-mesh.
%
% we need L to be defined as a sparse matrix,
% so we need arrays I and J to do:
%
%        L = sparse(I,J,v);
%
% where 'v' is special.
% ------------------------------------------------------------------------------
% we assume we already ran:
%       [n_ij,n_IJ] = nIJ(n_g2m,neigh_type);
% ------------------------------------------------------------------------------
% build I and J.
%
% both I and J are of size n_IJ by 1.
% I denotes the row entries, and J the column entries.
% for each node i:
% I needs to have n_ij(i)+1 consecutive entries with the number i.
% in the same place as these n_ij(i)+1 consecutive entries,
% J needs to have the number i in the first entry,
% and then each neighbor of i (in the graph) in the subsequent entries.
% ------------------------------------------------------------------------------
I = zeros(n_IJ,1);
J = zeros(n_IJ,1);
il = 1;
il_= 0;
for i_g2m=1:n_g2m
    % -- end of line
    il_ = il_ + n_ij(i_g2m)+1;
    % -- in line
    % I
    for ii=il:il_
        I(ii) = i_g2m;
    end
    % J
    J(il) = i_g2m;
    i_nei = 1;
    ii    = il+1;
    while (i_nei <= 6)
        if (neigh_graph(i_g2m,i_nei) ~= 0)
            J(ii) = neigh_graph(i_g2m,i_nei);
            ii = ii+1;
        end
        i_nei = i_nei + 1;
    end
    % -- begining of next line
    il = il_ + 1;
end
end
