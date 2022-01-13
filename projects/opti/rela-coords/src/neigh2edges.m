function [edges,nedges] = neigh2edges(np,max_neigh,neigh)
% diego domenzain
% ------------------------------------------------------------------------------
nedges = 0;
for ip=1:np
  ip_=1;
  while (ip_<=max_neigh && neigh(ip,ip_)~=0)
    if (ip<neigh(ip,ip_))
      nedges = nedges+1;
    end
    ip_ = ip_+1;
  end
end

edges = zeros(nedges,2,'uint32');
iedge = 1;
for ip=1:np
  ip_=1;
  while (ip_<=max_neigh && neigh(ip,ip_)~=0)
    if (ip<neigh(ip,ip_))
      edges(iedge,1) = ip;
      edges(iedge,2) = neigh(ip,ip_);
      iedge = iedge+1;
    end
    ip_ = ip_+1;
  end
end
end
