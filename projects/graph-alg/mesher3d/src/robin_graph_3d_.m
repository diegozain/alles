function robin_graph = robin_graph_3d_(n_g2m,nx,ny,nz,neigh_type)
% diego domenzain
% ? 2021
% ------------------------------------------------------------------------------
% neigh_type : row indexes are graph nodes.
%              row entries are the type of neighbor for that node.
%                   2  6
%                   | /
%              3 -- i -- 1
%                 / |
%                5  4
%
%              that is, columns 1, 2, 3, 4, 5 and 6,
%              represent neighbors right, up, left, down, front and back.
% ------------------------------------------------------------------------------
% get size of robin nodes
nprobin = 0;
for i_g2m=1:n_g2m
  i_ri = neigh_type(i_g2m,1);
  i_up = neigh_type(i_g2m,2);
  i_le = neigh_type(i_g2m,3);
  i_do = neigh_type(i_g2m,4);
  i_fw = neigh_type(i_g2m,5);
  i_bw = neigh_type(i_g2m,6);

  if (i_ri==-1)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
  if (i_up==-1)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
  if (i_le==-1)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
  if (i_do==-1)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
  if (i_fw==-1)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
  if (i_bw==-1)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
end
% ------------------------------------------------------------------------------
robin_mesh = zeros(nprobin,1,'uint32');

for i_g2m=1:nprobin
  i_ri = neigh_type(i_g2m,1);
  i_up = neigh_type(i_g2m,2);
  i_le = neigh_type(i_g2m,3);
  i_do = neigh_type(i_g2m,4);
  i_fw = neigh_type(i_g2m,5);
  i_bw = neigh_type(i_g2m,6);

  if (i_ri==-1)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin) = iyxz + nz;
  end
  if (i_up==-1)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin) = iyxz - 1;
  end
  if (i_le==-1)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin) = iyxz - nz;
  end
  if (i_do==-1)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin) = iyxz + 1;
  end
  if (i_fw==-1)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin) = iyxz - nz*nx;
  end
  if (i_bw==-1)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin) = iyxz + nz*nx;
  end
end
% ------------------------------------------------------------------------------
robin_mesh = unique(robin_mesh);
nprobin = numel(robin_mesh);
% ------------------------------------------------------------------------------
robin_mesh_xyz = zeros(nprobin,3,'uint32');
for iprobin=1:nprobin
  iyxz = robin_mesh(iprobin);
  [ix,iy,iz] = get_ixyz(iyxz,nx,ny,nz);
  robin_mesh_xyz(iprobin,:) = [ix,iy,iz];
end
end
