function neigh_graph = neigh_graph_(neigh_mesh,mesh2graph,n_g2m)
% diego domenzain
% April 2021 @ Colorado School of Mines
% 
% ------------------------------------------------------------------------------
% neigh_graph : row indexes are graph nodes.
%               row entries are neighbors of that node, in the graph.
% ------------------------------------------------------------------------------
neigh_graph = zeros(n_g2m,4);
for i_g2m = 1:n_g2m
    for i_nei = 1:4
        if (neigh_mesh(i_g2m,i_nei) ~= 0)
            neigh_graph(i_g2m,i_nei) = mesh2graph(neigh_mesh(i_g2m,i_nei));
        end
    end
end
end