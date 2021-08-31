function neigh_graph = neigh_graph_3d_(neigh_mesh,mesh2graph,n_g2m)
% diego domenzain
% August 2021
%
% ------------------------------------------------------------------------------
% neigh_graph : row indexes are graph nodes.
%               row entries are neighbors of that node, in the graph.
% ------------------------------------------------------------------------------
neigh_graph = zeros(n_g2m,6,'uint32');
for i_g2m = 1:n_g2m
    for i_nei = 1:6
        if (neigh_mesh(i_g2m,i_nei) ~= 0)
            neigh_graph(i_g2m,i_nei) = mesh2graph(neigh_mesh(i_g2m,i_nei));
        end
    end
end
end
