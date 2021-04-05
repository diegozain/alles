function [n_ij,n_IJ] = nIJ(n_g2m,neigh_type)
% diego domenzain
% April 2021 @ Colorado School of Mines
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
n_ij = zeros(n_g2m,1);
n_IJ = 0;

for i_g2m=1:n_g2m
    n_i = 0;
    for i_nei=1:4
        if (neigh_type(i_g2m,i_nei) > 0)
            n_i = n_i + neigh_type(i_g2m,i_nei);
        end
    end
    n_ij(i_g2m) = n_i;
    n_IJ = n_IJ + (n_i + 1);
end

end