function [I,J] = IJ_graph(ngraph,nei_max,n_ij,n_IJ,neigh_graph)
% diego domenzain
% december 2021 @ Tepoztlan
%
% ------------------------------------------------------------------------------
% let L be laplacian matrix of the 🍇.
%
% we need L to be defined as a sparse matrix,
% so we need arrays I and J to do:
%
%        L = sparse(I,J,v);
%
% where 'v' is special.
% ------------------------------------------------------------------------------
% we assume we already ran:
%       [n_ij,n_IJ] = nIJ_graph(ngraph,nei_max,neigh_graph);
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
for ig=1:ngraph
    % -- end of line
    il_ = il_ + n_ij(ig)+1;
    % -- in line
    % I
    for ii=il:il_
        I(ii) = ig;
    end
    % J
    J(il) = ig;
    i_nei = 1;
    ii    = il+1;
    while (i_nei<=nei_max)
        if (neigh_graph(ig,i_nei) ~= 0)
            J(ii) = neigh_graph(ig,i_nei);
            ii = ii+1;
        end
        i_nei = i_nei + 1;
    end
    % -- begining of next line
    il = il_ + 1;
end
end
