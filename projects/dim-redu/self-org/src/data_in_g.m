function d_in_g = data_in_g(g,d,w)
% which data points are near which graph points?
n_nodes = size(g,1);
nd = size(d,2);
d_in_g = zeros(nd,1);
for id=1:nd % parforable
  d_=d(:,id);
  dist_ = zeros(n_nodes,1);
  for ig=1:n_nodes % parforable
    dist_(ig) = norm(d_-w(:,ig));
  end
  [~,mini] = min(dist_);
  d_in_g(id) = mini;
end
end