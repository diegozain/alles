function [w,u_mat,d_in_g] = selforgmapi(d,g,w,n_iter)
% diego domenzain
% Boise State University
% ---------------------------------------------------------------------------
% builds self-organizing map from:
% data points d and graph g.
% ---------
% g -> is an incidence relation, e.g. (but not limited to):
%		g = graph_grid(ny,nx);
%		that is, a cube (ny x nx x 5) where:
%		the first (ny x nx) plane labels of the nodes,
%		the second (ny x nx) plane labels the right neighbors,
%		etc.
%		g has to be of size (#of-nodes by max-#of-neighbors)
% d -> is data where rows are attributes and 
% 	columns are data points.
% w -> are initial weights on graph g, e.g. geometrical coordinates in R^n. 
%		Rows are entries of R^n, columns are nodes in the grid.
%		For example, 
%		w = w_in_grid(v,u,ny,nx);
%		where v and u are the vectors spanning a plane in R^n.
% ----------
% u_mat -> is a u-matrix: a matrix in form of a list 
%		of size (#of-nodes by 1) where rows are nodes indexed 
%		by their number and entry of row i is the average distance
%		of node #i to its neighbors. For example,
%		u_mat = u_matrix(g,w);
% d_in_g -> is a matrix of size (#of-data-pts by 1) 
%		where each entry answers:
%		"which data points are near which graph points?", that is
%		entry of row i is the graph node closest to data point i (euclid).
%	  For example,
%		d_in_g = data_in_g(g,d,w);
% ------------------------------------------------------------------------------
% data points are columns, attributes are rows,
% like a boring excel sheet transposed.
[n_atributes,nd] = size(d);
n_nodes = size(g,1);
% -----------
std_ = 8;
amp = 2;
tau = 1;
% ------------------------------------------------------------------------------
%
% main loop
%
% ------------------------------------------------------------------------------
iter = 0;
choose_list = randperm(nd);
while iter<n_iter
  % choose data point
  if numel(choose_list)==0
    choose_list = randperm(nd);
  end
  d_ = d(:,choose_list(1));
  choose_list(1) = [];
  % initialize "best matching unit" (bmu)
  bmu = zeros(n_nodes,1);
  for i_=1:n_nodes % parforable
    bmu(i_) = norm(d_-w(:,i_));
  end
  % get bmu: bmu has to be w(:,bmu)
  [~,bmu] = min(bmu);
  % get neighbors of bmu,
  % non-zero entries of g at node bmu
  nei = find(g(bmu,:));
  % these are the actual neighboors
  nei = g(bmu,nei);
  % AND bmu itself
  nei = [nei bmu];
  n_nei = numel(nei);
  % update neighbors
  for i_=1:n_nei % parforable
    % w(:,nei(i_)) = w(:,nei(i_)) +...
    %  thet(w(:,bmu),w(:,nei(i_)),std_,iter)*alp(iter,amp,tau)*(d_-w(:,nei(i_)));
    % r=rand(n_atributes,1);
    % w(:,nei(i_)) = w(:,nei(i_)) + 0.4*(d_-w(:,nei(i_))) + 0.01*(2*(r-mean(r)));
    w(:,nei(i_)) = w(:,nei(i_)) + 0.4*(d_-w(:,nei(i_)));
  end
  iter = iter+1;
end
% --------------
% now that the main loop is done,
% we have all nodes in the graph covering the data points in data space:
% node i is in position w(:,i).
% ---------------
% now get cute bw plot of distances between neighboring weights,
% aka u-matrix.
u_mat = u_matrix(g,w);
% ---------------
% which data points are near which graph points?
d_in_g = data_in_g(g,d,w);
end
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
function t = thet(u,v,std_,iter)
  dist = norm(u-v);
  t = exp(-dist/(2*std_^2) );
end
function a = alp(iter,amp,tau)
  a = amp*exp(- iter );
end