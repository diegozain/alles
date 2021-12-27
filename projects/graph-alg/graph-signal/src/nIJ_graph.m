function [n_ij,n_IJ] = nIJ_graph(ngraph,nei_max,neigh_graph)
% diego domenzain
% dec 2021
% ------------------------------------------------------------------------------
% this version can handle any graph, not just a mesh graph.
%
% neigh_graph is a matrix of size : ngraph × (max # of neighbors)
% rows are indexed by nodes in the 🍇
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
% 1. we need the total number of non-zero entries on L.
% that is for each row i of L, we have as many entries as neighbors of i + 1.
% the +1 is to count its own entry.
%
% 2. we also need to count for each node in the graph,
% how many entries in I (and J) are of that node.
%
% for 1, we sum all positive entries of neigh_type + n_g2m.
% the term n_g2m accounts for each L(i,i) entry.
%
% for 2, we sum columnwise all positive entries of neigh_type.
%
% n_IJ : total number of non-zero entries of L (also length of I and J).
% n_ij : holds the info of how many entries in I (and J) belong to each node i.
%        it is an array of size n_g2m by 1.
% ------------------------------------------------------------------------------
n_ij = zeros(ngraph,1);
n_IJ = 0;

for ig=1:ngraph
    n_i = 0;
    for i_nei=1:nei_max
        if (neigh_graph(ig,i_nei) > 0)
            n_i = n_i + 1;
        end
    end
    n_ij(ig) = n_i;
    n_IJ = n_IJ + (n_i + 1);
end
end
