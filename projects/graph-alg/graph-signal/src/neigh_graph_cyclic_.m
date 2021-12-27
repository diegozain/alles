function neigh_graph = neigh_graph_cyclic_(ngraph)
% diego domenzain
% dec 2021
% ------------------------------------------------------------------------------
% all cyclic graphs are 2-regular
neigh_graph = zeros(ngraph,2,'uint32');

% now fill in 'neigh_graph'
for ig=1:ngraph
  nei_ = mod(ig+1,ngraph);
  if nei_ == 0
    nei_ = ngraph;
  end
  neigh_graph(ig,1) = nei_;

  nei_ = mod(ig-1,ngraph);
  if nei_ == 0
    nei_ = ngraph;
  end
  neigh_graph(ig,2) = nei_;
end
end
