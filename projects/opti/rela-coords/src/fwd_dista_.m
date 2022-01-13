function dista = fwd_dista_(np,nedges,xy,edges)
% diego domenzain
dista = zeros(nedges,1);
for iedge=1:nedges
  dista(iedge) = ( xy(edges(iedge,1)) - xy(edges(iedge,2)) )^2 + ( xy(edges(iedge,1)+np) - xy(edges(iedge,2)+np) )^2;
end
end
