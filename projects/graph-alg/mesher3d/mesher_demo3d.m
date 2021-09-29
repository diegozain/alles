clear
close all
clc
clc
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
%                            mesh generation
%
% given a rectangular gridded mesh where only a certain region is wanted for
% modeling, how do we extract this wanted region?
%
% Answer:
%
% 1. find only the nodes that matter in the domain,
% 2. find their neighbors,
% 3. find what type of neighbors they are.
%
% The words **mesh** or **mesh-grid** refer to the existing rectangular gridded mesh.
%
% The word **graph** refers to the region of interest inside the mesh-grid.
%
% Expanded answer:
%
% 1. The only nodes that matter in the domain are given in these two column vectors,
%
% graph2mesh : indexes are graph nodes, entries are mesh nodes
% mesh2graph : indexes are mesh nodes, entries are graph nodes
%
% graph2mesh:
%
% | |
% | |
% | | size = # of nodes that matter by 1
% | |
% | |
% | |
%
% mesh2graph:
%
% | |
% | |
% | | size = # of nodes in the mesh by 1
% | |
% | |
% | |
% | |
% | |
%
% 2. Their neighbors are given by these two matrices,
%
% neigh_mesh  : row indexes are graph nodes. Row entries are neighbors of that node, in the mesh.
% neigh_graph : row indexes are graph nodes. Row entries are neighbors of that node, in the graph.
%
% 3. Their neighbor type is given by this matrix,
%
% neigh_type : row indexes are graph nodes. Row entries are the type of neighbor for that node.
%
% neigh_mesh, neigh_graph, and neigh_type are all of size:
%
% |        |
% |        |
% |        | size = # of nodes that matter by max # of neighbors in the mesh
% |        |
% |        |
% |        |
%
% since we are assuming a 2D mesh-grid, max # of neighbors in the mesh is 4.
%
% the way the neigh_mesh, neigh_graph, and neigh_type refer to the neighbors of node 'i' (i.e. in row 'i') is by,
%
%      2  6
%      | /
% 3 -- i -- 1
%    / |
%   5  4
%
% that is, columns 1, 2, 3, 4, 5 and 6,
% represent neighbors right, up, left, down, front and back.
% ------------------------------------------------------------------------------
% -- setup a simple example
ay1 = [0 0 0 1 1 0 0; 0 1 1 1 1 0 0; 1 1 1 1 1 1 1; 1 1 0 0 0 1 1; 1 1 1 1 1 1 1];
ay2 = [0 0 0 1 1 0 0; 0 1 1 1 1 0 0; 1 1 1 1 1 1 1; 1 1 0 0 0 1 1; 1 1 1 1 1 1 1];
ay3 = [0 0 0 1 1 0 0; 0 1 1 1 1 0 0; 1 1 1 1 1 1 1; 1 1 0 0 0 1 1; 1 1 1 1 1 1 1];
ay4 = [0 0 0 1 1 0 0; 0 1 1 1 1 0 0; 1 1 1 1 1 1 1; 1 1 0 0 0 1 1; 1 1 1 1 1 1 1];
a_ = cat(3,ay1,ay2,ay3,ay4);
a_ = permute(a_,[3,2,1]);

[ny,nx,nz]=size(a_);
% ------------------------------------------------------------------------------
% this is only for easy reference:
a_index = 1:(nx*ny*nz);
a_index = reshape(a_index,[nz,nx,ny]);
a_index = permute(a_index,[3,2,1]);
% this is the mask 😷 in the mesh
amat3d_mask = zeros(ny*nx*nz,4,'uint32');
% this is the 3d mesh 🎲 with numering scheme
amat3d_mesh= zeros(ny*nx*nz,4,'uint32');

% this is the actual 3d matrix 🎲 in the mesh in just one column
a = zeros(ny*nx*nz,1);
for iyxz = 1:ny*nx*nz
  % get x,y,z coordinate
  [ix,iy,iz] = get_ixyz(iyxz,nx,ny,nz);

  % this is just for vis 🎨
  amat3d_mask(iyxz,:) = [iy,ix,iz,a_(iy,ix,iz)];
  amat3d_mesh(iyxz,:) = [iy,ix,iz,a_index(iy,ix,iz)];

  % 3d matrix 🎲 in the mesh in just one column
  a(iyxz,:) = a_(iy,ix,iz);
end
% ------------------------------------------------------------------------------
% get the number of nodes in the graph
n_g2m = n_g2m_3d_(a,nx,ny,nz);
% ------------------------------------------------------------------------------
% make two dictionaries,
% graph2mesh : indexes are graph nodes, entries are mesh nodes
% mesh2graph : indexes are mesh nodes, entries are graph nodes
[graph2mesh,mesh2graph] = g2m_m2g_3d(a,nx,ny,nz,n_g2m);
% ------------------------------------------------------------------------------
% each node has neighbors.
% neigh_mesh : row indexes are graph nodes.
%              row entries are neighbors of that node, in the mesh.
neigh_mesh = neigh_mesh_3d_(a,nx,ny,nz,n_g2m,graph2mesh);
% ------------------------------------------------------------------------------
% each node has neighbors.
% neigh_graph : row indexes are graph nodes.
%               row entries are neighbors of that node, in the graph.
neigh_graph = neigh_graph_3d_(neigh_mesh,mesh2graph,n_g2m);
% ------------------------------------------------------------------------------
% each node has a special type in a mesh-grid.
% neighbors that are (in the mesh):
% zero, non-zero, and next to the limits of the mesh.
%
% In the case the mesh is a slice of the earth, and we are doing this for a PDE,
% this translates to:
%
% neighbors that are non-zero = the pde (inner-node)
% neighbors that are zero = neumann bc
% neighbors that are next to the limits of the mesh = robin bc
%
% neigh_type : row indexes are graph nodes.
%              row entries are the type of neighbor for that node.
%
% we define : (type,BC) = (1,inner) (-1,neumann) (0,robin)
% ------------------------------------------------------------------------------
neigh_type = neigh_type_3d_(a,nx,ny,nz,n_g2m,graph2mesh);
% ------------------------------------------------------------------------------
%
%                              🎨 vis 🎨
%
% ------------------------------------------------------------------------------
% get robin nodes
[robin_xyz,robin_mesh,robin_mesh_] = robins_3d_(n_g2m,nx,ny,nz,graph2mesh,neigh_type);
% ------------------------------------------------------------------------------
nprobin_ = size(robin_mesh_,1);
% translate to coordinates in the mesh 🎲
robin_xyz_ = zeros(nprobin_,4,'uint32');
for iprobin_=1:nprobin_
  iyxz = robin_mesh_(iprobin_,1);
  [ix,iy,iz] = get_ixyz(iyxz,nx,ny,nz);
  robin_xyz_(iprobin_,1:3) = [iy,ix,iz];
  robin_xyz_(iprobin_,4) = robin_mesh_(iprobin_,2);
end
% ------------------------------------------------------------------------------
% this is the 3d graph 🍇 with indexing scheme
amat3d_graph = zeros(n_g2m,4,'uint32');
for i_g2m=1:n_g2m
  % get x,y,z coordinate
  iyxz = graph2mesh(i_g2m);
  [ix,iy,iz] = get_ixyz(iyxz,nx,ny,nz);
  % 3d graph 🍇 with indexing scheme
  amat3d_graph(i_g2m,:) = [iy,ix,iz,i_g2m];
end
% ------------------------------------------------------------------------------
figure;
scatter3(amat3d_mask(:,2),amat3d_mask(:,1),amat3d_mask(:,3),500*abs(amat3d_mask(:,4))+1,amat3d_mask(:,4),'filled')
colormap('lines');
set(gca,'ZDir','reverse');
axis image;
axis tight;
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh-mask 😷')
simple_figure()

figure;
subplot(1,2,1);
scatter3(amat3d_mesh(:,2),amat3d_mesh(:,1),amat3d_mesh(:,3),500*ones(ny*nx*nz,1),amat3d_mesh(:,4),'filled')
colormap(rainbow2_cb(1));
hcb = colorbar('southoutside');
hcb.TickLength = 0;
xlabel(hcb,'Index in mesh-cube');
set(gca,'ZDir','reverse');
axis image;
axis tight;
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh 🎲')
simple_figure()

subplot(1,2,2);
scatter3(amat3d_graph(:,2),amat3d_graph(:,1),amat3d_graph(:,3),500*ones(n_g2m,1),amat3d_graph(:,4),'filled')
colormap(rainbow2_cb(1));
hcb = colorbar('southoutside');
hcb.TickLength = 0;
xlabel(hcb,'Index in graph');
set(gca,'ZDir','reverse');
axis image;
axis tight;
xlabel('x')
ylabel('y')
zlabel('z')
title('Graph 🍇')
simple_figure()

figure;
subplot(1,2,1);
scatter3(robin_xyz_(:,2),robin_xyz_(:,1),robin_xyz_(:,3),500*ones(nprobin_,1),robin_xyz_(:,4),'filled')
colormap(rainbow2_cb(1));
hcb = colorbar('southoutside');
hcb.TickLength = 0;
xlabel(hcb,'# of 🐦 entries');
set(gca,'ZDir','reverse');
axis image;
axis tight;
xlabel('x')
ylabel('y')
zlabel('z')
title('Robin in mesh 🎲')
simple_figure()
% ------------------------------------------------------------------------------
%
%     🔷 this part is for building the PDE operator L acting on the graph 🔺
%
% ------------------------------------------------------------------------------
% we need L to be defined as a sparse matrix,
% so we need arrays I and J to do:
%
%        L = sparse(I,J,v);
%
% where 'v' is special.
% ------------------------------------------------------------------------------
% 1. we need the total number of non-zero entries on L.
% that is for each row i of L, we have as many entries as neighbors of i + 1.
% the +1 is to count its own entry.
%
% 2. we also need to count for each node in the graph,
% how many entries in I (and J) are of that node.
%
% for 1, we sum all positive entries of neigh_type +n_g2m.
% the term +n_g2m accounts for all L(i,i) entries.
%
% for 2, we sum columnwise all positive entries of neigh_type.
%
% n_IJ : total number of non-zero entries of L (also length of I and J).
% n_ij : holds the info of how many entries in I (and J) belong to each node i.
%        it is an array of size n_g2m by 1.
% ------------------------------------------------------------------------------
[n_ij,n_IJ] = nIJ_3d(n_g2m,neigh_type);
% ------------------------------------------------------------------------------
% now we need to build I and J.
%
% both I and J are of size n_IJ by 1.
% I denotes the row entries, and J the column entries.
% for each node i:
% I needs to have n_ij(i)+1 consecutive entries with the number i.
% in the same place as these n_ij(i)+1 consecutive entries,
% J needs to have the number i in the first entry,
% and then each neighbor of i (in the graph) in the subsequent entries.
% ------------------------------------------------------------------------------
[I,J] = IJ_3d_(n_g2m,n_ij,n_IJ,neigh_graph);
% ------------------------------------------------------------------------------
%
%                         🐦 α's for robin bc 🐦
%
% ------------------------------------------------------------------------------
dx=0.5; % m
dy=0.5; % m
dz=0.5; % m

x=(0:(nx-1))*dx; x=x.';
y=(0:(ny-1))*dy; y=y.';
z=(0:(nz-1))*dz; z=z.';
% ------------------------------------------------------------------------------
%
% srcs_xyz  : (nsources) × (xyz) × (±) . indexes in the mesh cube 🎲
srcs_xyz = zeros(1,3,2,'uint32');
srcs_xyz(1,1,1) = 3; % x index of a
srcs_xyz(1,2,1) = 2; % y index of a
srcs_xyz(1,3,1) = 2; % z index of a

srcs_xyz(1,1,2) = 6; % x index of b
srcs_xyz(1,2,2) = 2; % y index of b
srcs_xyz(1,3,2) = 3; % z index of b
% ------------------------------------------------------------------------------
iy = srcs_xyz(1,2,1);
ix = srcs_xyz(1,1,1);
iz = srcs_xyz(1,3,1);
iyxz = (iy-1)*nx*nz + (ix-1)*nz + iz;
a_in_g = mesh2graph(iyxz);

iy = srcs_xyz(1,2,2);
ix = srcs_xyz(1,1,2);
iz = srcs_xyz(1,3,2);
iyxz = (iy-1)*nx*nz + (ix-1)*nz + iz;
b_in_g = mesh2graph(iyxz);
% ------------------------------------------------------------------------------
% 🐦
[robin_graph,robin_xyz,n_ar] = robins_3d(n_g2m,nx,ny,nz,graph2mesh,mesh2graph,neigh_type);
% α
alphas = get_alphas(x,y,z,srcs_xyz,robin_xyz);
% ------------------------------------------------------------------------------
nprobin= size(robin_xyz,1);
nsource= size(srcs_xyz,1);

subplot(1,2,2);
scatter3(robin_xyz(:,1),robin_xyz(:,2),robin_xyz(:,3),300*ones(nprobin,1),alphas,'filled')
colormap(rainbow2_cb(1));
hcb = colorbar('southoutside');
hcb.TickLength = 0;
xlabel(hcb,'α');
set(gca,'ZDir','reverse');
axis image;
axis tight;
xlabel('x')
ylabel('y')
zlabel('z')
title('α 🐦')
simple_figure()
% ------------------------------------------------------------------------------
%
%
%                            🔷  build L  🔺
%
%
% ------------------------------------------------------------------------------
% declare σ
sig=ones(n_g2m,1);
% ------------------------------------------------------------------------------
V = dcipL3d(n_g2m,n_ij,n_IJ,I,J,neigh_mesh,graph2mesh,robin_graph,alphas,n_ar,sig,x,y,z);
% ------------------------------------------------------------------------------
L = sparse(I,J,V);
% ------------------------------------------------------------------------------
% this is how we would solve for the electric potential 'u'
%       -∇⋅σ ∇ u = s
% ------------------------------------------------------------------------------
% source
s=zeros(n_g2m,1);
s(a_in_g)= 1;
s(b_in_g)=-1;
% solve
u=L\s;
% ------------------------------------------------------------------------------
%
%                               🎨 vis 🎨
%
% ------------------------------------------------------------------------------
% this is in the 3d graph 🍇
u3d = zeros(n_g2m,4);
for i_g2m=1:n_g2m
  % get x,y,z coordinate
  iyxz = graph2mesh(i_g2m);
  [ix,iy,iz] = get_ixyz(iyxz,nx,ny,nz);
  % 🎨
  u3d(i_g2m,1:3)= [iy,ix,iz];
  u3d(i_g2m,4)  = u(i_g2m);
end

figure;
% scatter3(u3d(:,2),u3d(:,1),u3d(:,3),500*ones(n_g2m,1),u3d(:,4),'filled')
scatter3(u3d(:,2),u3d(:,1),u3d(:,3),500*abs(u3d(:,4)),u3d(:,4),'filled')
colormap(rainbow2_cb(1));
hcb = colorbar('southoutside');
hcb.TickLength = 0;
xlabel(hcb,'(🌝)');
set(gca,'ZDir','reverse');
axis image;
axis tight;
xlabel('x')
ylabel('y')
zlabel('z')
title('ϕ 🔌')
simple_figure()
% ------------------------------------------------------------------------------
