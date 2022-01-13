function Jac = jac_dista_(np,nedges,max_neigh,neigh,xy)
% diego domenzain
% ------------------------------------------------------------------------------
% Jac (2⋅(# of points) × # of edges)
%
% because edges not spanning from p(ip) are zero under ∂x_ip & ∂y_ip,
%
% Jac( ip,nnz(neigh(ip,:)) ) = [ ∂x_ip (distances of edges spanning from ip) ]
% Jac( ip+np,nnz(neigh(ip,:)))=[ ∂y_ip (distances of edges spanning from ip) ]
%
% ------------------------------------------------------------------------------
Jac = zeros(2*np,nedges);
iedge = 1;
for ip = 1:np
  ip_=1;
  while (ip_<=max_neigh && neigh(ip,ip_)~=0)
    if (ip<neigh(ip,ip_))
      % x
      Jac(ip,iedge)   = 2*( xy(ip) - xy(neigh(ip,ip_)) );
      % y
      Jac(ip+np,iedge)= 2*( xy(ip+np) - xy(neigh(ip,ip_)+np) );

      % neigh of p(ip)
      Jac(neigh(ip,ip_),iedge)   = 2*( xy(neigh(ip,ip_)) - xy(ip) );
      Jac(neigh(ip,ip_)+np,iedge)= 2*( xy(neigh(ip,ip_)+np) - xy(ip+np));

      iedge = iedge+1;
    end
    ip_ = ip_+1;
  end
end
end
