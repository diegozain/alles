function [u,P] = fast_marcher(A,c,s,h)
% given velocity c among nodes,
% an adjacency matrix A,
% and source position s,
% return time field u, 
% and **^^-- magic --^^** vector P.

% get size of field
[nx,nz] = size(c);
c = c(:);
nu = nx*nz;
% list of nodes:
%
% [ not-accepted ]
Q = sparse(nu,1);
% [ path ]
P = sparse(nu,1);

% initialize u 
u = inf(nu,1);
% set field on source to zero
u(s) = 0;
Q(s,1) = 1;
% initialize for all nodes
for i_=1:nu
  % skip sources
  if u(i_) == 0
    continue;
  end
  % loop over neighbors of i_ and
  % update field u
  i_neighbors = find( A(:,i_) );
  [u,~] = u_update( i_,i_neighbors,u,c,h,nx );
end

% loop over not-accepted nodes
while sum( Q(:,1) ) < nu

  % find node with least value u from not-accepted
  % and accept
  q = find( ~Q(:,1) );
  [~,i_] = min( u( q ) );
  i_ = q( i_ );
  Q(i_,1) = 1;

  % loop over neighbors of i_ and
  % update field u
  i_neighbors = find( A(:,i_) );
  for j__ = 1:numel( i_neighbors )
    j_ = i_neighbors( j__ );
    j_neighbors = find( A(:,j_) );
    [u,u__] = u_update( j_,j_neighbors,u,c,h,nx );
    if u(j_) < u__
      P(j_) = i_;
    end
  end
end
end

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

function [u,u__] = u_update(i_,i_neighbors,u,c,h,nx)

j_x = [i_neighbors(i_neighbors==i_-1); i_neighbors(i_neighbors==i_+1)];
j_z = [i_neighbors(i_neighbors==i_-nx); i_neighbors(i_neighbors==i_+nx)];

ux = min( u(j_x) ); 
uz = min( u(j_z) );


  % eikonal update
  if abs(ux-uz) > (h/c(i_))
    u_ = min(ux,uz) + h/c(i_);
  else
    u_ = 0.5*( (ux+uz) + sqrt( (ux+uz)^2 - 2*(ux^2 + uz^2 - (h/c(i_))^2)  ) );
  end
  
  u__ = u(i_);
  if u_ < u__
    u(i_) = u_;
  end

end
