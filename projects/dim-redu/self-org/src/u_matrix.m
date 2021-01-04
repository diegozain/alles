function u_mat = u_matrix(g,w)
% get u-matrix from graph g and weights w
n_nodes = size(g,1);
u_mat = zeros(n_nodes,1);
for i_=1:n_nodes % parforable
  % get neighbors of node i_
  nei = find(g(i_,:));
  % these are the actual neighboors
  nei = g(i_,nei);
  n_nei = numel(nei);
  dist_i = zeros(n_nei,1);
  for j_=1:n_nei % parforable
    dist_i(j_) = norm(w(:,i_)-w(:,nei(j_)));
  end
  % get average distance
  u_mat(i_) = sum(dist_i)/n_nei;
end
u_mat = u_mat / max(u_mat(:));
end