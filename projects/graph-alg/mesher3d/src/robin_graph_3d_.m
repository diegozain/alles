function [robin_xyz,robin_mesh,robin_graph] = robin_graph_3d_(n_g2m,nx,ny,nz,graph2mesh,mesh2graph,neigh_type)
% diego domenzain
% ? 2021
% ------------------------------------------------------------------------------
% ⚫ get nodes that are of type Robin in the mesh.
%
% when 'neigh_type' is built, it is initialized as all robin and then corrects
% for each node which type it is.
% this is because getting the robin nodes themselves is tricky.
%
% so here we need to count robin nodes inside 'neigh_type', but that itself is
% tricky because some robin nodes can appear more than once in this list.
%
% 1. count how many nodes are of type robin in 'neigh_type':
%                                                        ▶️ nprobin.
% 2. make a list of size nprobin and fill it in with the actual node id.
% 3. each node in this list carries info on which neighbors are robin:
%                                                        ▶️ robin_mesh.
% 4. take this unique list and also save their xyz index coordinates:
%                                                        ▶️ robin_xyz.
% ------------------------------------------------------------------------------
% ▶️
% robin_mesh : robin nodes in the mesh 🎲. of size nprobin × 2.
%              in the second column you find which side is robin.
%              for example, a corner node will be repeated 3 times in robin_mesh
%              and the second column will perhaps read 3,4,5, meaning
%              left, down, front neighbors are robin.
% robin_xyz  : robin nodes in the mesh cube 🎲. of size nprobin × 4.
%              the 4rth column is the second column of robin_mesh.
% robin_graph: robin nodes in the graph 🍇. of size nprobin × 1.
% ------------------------------------------------------------------------------
% ◀️
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

  if (i_ri==0)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
  if (i_up==0)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
  if (i_le==0)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
  if (i_do==0)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
  if (i_fw==0)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
  if (i_bw==0)
    iyxz = graph2mesh(i_g2m);
    nprobin = nprobin + 1;
  end
end
% ------------------------------------------------------------------------------
robin_mesh = zeros(nprobin,2,'uint32');
iprobin=1;
for i_g2m=1:n_g2m
  i_ri = neigh_type(i_g2m,1);
  i_up = neigh_type(i_g2m,2);
  i_le = neigh_type(i_g2m,3);
  i_do = neigh_type(i_g2m,4);
  i_fw = neigh_type(i_g2m,5);
  i_bw = neigh_type(i_g2m,6);

  if (i_ri==0)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin,1) = iyxz;
    robin_mesh(iprobin,2) = 1;
    iprobin = iprobin + 1;
  end
  if (i_up==0)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin,1) = iyxz;
    robin_mesh(iprobin,2) = 2;
    iprobin = iprobin + 1;
  end
  if (i_le==0)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin,1) = iyxz;
    robin_mesh(iprobin,2) = 3;
    iprobin = iprobin + 1;
  end
  if (i_do==0)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin,1) = iyxz;
    robin_mesh(iprobin,2) = 4;
    iprobin = iprobin + 1;
  end
  if (i_fw==0)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin,1) = iyxz;
    robin_mesh(iprobin,2) = 5;
    iprobin = iprobin + 1;
  end
  if (i_bw==0)
    iyxz = graph2mesh(i_g2m);
    robin_mesh(iprobin,1) = iyxz;
    robin_mesh(iprobin,2) = 6;
    iprobin = iprobin + 1;
  end
end
% ------------------------------------------------------------------------------
% translate to the graph 🍇
robin_graph = mesh2graph(robin_mesh(:,1));
% ------------------------------------------------------------------------------
% translate to coordinates in the mesh 🎲
robin_xyz = zeros(nprobin,4,'uint32');
for iprobin=1:nprobin
  iyxz = robin_mesh(iprobin,1);
  [ix,iy,iz] = get_ixyz(iyxz,nx,ny,nz);
  robin_xyz(iprobin,1:3) = [iy,ix,iz];
  robin_xyz(iprobin,4) = robin_mesh(iprobin,2);
end
end
