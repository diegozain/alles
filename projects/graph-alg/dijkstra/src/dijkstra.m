function [u,P] = dijkstra(c,s)
% given velocity c among nodes 
% (c is a weighted adjacency matrix) and
% source position s,
% return time field u,
% and **^^magic^^** vector P.

% get size of field
[nu,~] = size(c);

% initialize u 
u = inf(nu,1);
% set field on source to zero
u(s) = 0;

% list of nodes:
%
% [ not-accepted ]
Q = sparse(nu,1);
% [ path ]
P = sparse(nu,1);

% loop over not-accepted nodes
while sum( Q(:,1) ) < nu
  
  % find node with least value u from not-accepted
  % and accpet
  q = find( ~Q(:,1) );
  [~,i_] = min( u( q ) );
  i_ = q( i_ );
  Q(i_,1) = 1;

  % loop over neighbors of i_ and
  % update field u
  i_neighbors = find( c(:,i_) );
  for j__ = 1:numel( i_neighbors )
    j_ = i_neighbors( j__ );
    u_ = u_update( i_,j_,u,c );
    if u_ < u(j_)
      u(j_) = u_;
      P(j_) = i_;
    end
  end
end
end

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

function u_ = u_update(i_,j_,u,c)
% simple distance field
u_ = u(i_) + c(i_,j_);
end