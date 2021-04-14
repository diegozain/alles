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
%      2
%      |
% 3 -- i -- 1
%      |
%      4
% 
% that is, columns 1, 2, 3, and 4 represent neighbors right, up, left and down.
% ------------------------------------------------------------------------------
% setup a simple example
% a = [0 0 0 1 1 0 0; 0 1 1 1 1 0 0; 1 1 1 1 1 1 1; 1 1 0 0 0 1 1];
a = [0 0 0 1 1 0 0; 0 1 1 1 1 0 0; 1 1 1 1 1 1 1; 1 1 0 0 0 1 1; 1 1 1 1 1 1 1];
[nz,nx]=size(a);
% ------------------------------------------------------------------------------
% this is only for easy reference:
a_index = 1:(nx*nz);
a_index = reshape(a_index,[nz,nx]);
% ------------------------------------------------------------------------------
% get the number of nodes in the graph 
n_g2m = n_g2m_(a,nx,nz);
% ------------------------------------------------------------------------------
% make two dictionaries,
% graph2mesh : indexes are graph nodes, entries are mesh nodes
% mesh2graph : indexes are mesh nodes, entries are graph nodes
[graph2mesh,mesh2graph] = g2m_m2g(a,nx,nz,n_g2m);
% ------------------------------------------------------------------------------
% each node has neighbors.
% neigh_mesh : row indexes are graph nodes.
%              row entries are neighbors of that node, in the mesh. 
neigh_mesh = neigh_mesh_(a,nx,nz,n_g2m,graph2mesh);
% ------------------------------------------------------------------------------
% each node has neighbors.
% neigh_graph : row indexes are graph nodes.
%               row entries are neighbors of that node, in the graph.
neigh_graph = neigh_graph_(neigh_mesh,mesh2graph,n_g2m);
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
neigh_type = neigh_type_(a,nx,nz,n_g2m,graph2mesh);
% ------------------------------------------------------------------------------
% 
%       this part is for building the PDE operator L acting on the graph
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
[n_ij,n_IJ] = nIJ(n_g2m,neigh_type);
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
[I,J] = IJ_(n_g2m,n_ij,n_IJ,neigh_graph);
% ------------------------------------------------------------------------------
V = zeros(n_IJ,1);
il = 1;
il_= 0;
for i_g2m=1:n_g2m
    % -- end of line
    il_ = il_ + n_ij(i_g2m)+1;
    % -- in line
    % V
    ith = 0;
    for ii=(il+1):il_
        % sum for ith entry
        ith = ith +1;
        % neighbors entries
        V(ii) = -1;
        % NOTE: the +1 and -1 should be replaced (each) with a product of:
        % 1. material properties (i.e. harmonic averages of conductivity),
        % 2. and dx_j / dz_j terms.
        % 
        % to do so, use J.
        % for example, inside this loop:
        % J(il) gives the ith node,
        % J(ii) gives neighbor of ith node.
        % 
        % the harmonic average would be:
        % sig_ij_ = ( 2*sig(J(il)) * sig(J(ii)) ) / ( sig(J(il)) + sig(J(ii)) )
    end
    % robin nodes are summed to 'ith'
    for i_nei=1:4
        if (neigh_type(J(il),i_nei) == 0)
            ith = ith +1;
            % NOTE: the +1 should be replaced by the appropriate entry.
            % J(il) gives the ith node,
            % neigh_type(J(il),i_nei) gives neighbor of ith node.
            % 
            % WARNING: be careful with the corners!!
            % 'bottom' corners are characterized by having two 0 entries in:
            %     neigh_type(J(il),1:4)
            % 'top' corners are characterized by having one 0 and one -1 in:
            %     neigh_type(J(il),1:4)
        end
    end
    % ith entry
    V(il) = ith;
    % -- begining of next line
    il = il_ + 1;
end
% ------------------------------------------------------------------------------
L = sparse(I,J,V);
% ------------------------------------------------------------------------------
% this is how we would solve for the electric potential 'u'
%       -∇⋅σ ∇ u = s
% ------------------------------------------------------------------------------
% source
s=zeros(n_g2m,1);
s(4)=1;
s(19)=-1;
% solve
u=L\s;
% for plotting
u_=nan(nz,nx);
u_(graph2mesh)=u;
u_plot=[u_;nan(1,nx)];
u_plot=[u_plot,nan(nz+1,1)];
% ------------------------------------------------------------------------------
% 
%                              visualize results
% 
% ------------------------------------------------------------------------------
fprintf(' -------------- mesh under consideration -------------- \n')
a

fprintf(' -------------- indexes of mesh ----------------------- \n')
a_index

fprintf(' -------------- graph -> mesh ------------------------- \n')
graph2mesh

fprintf(' -------------- mesh -> graph ------------------------- \n')
mesh2graph

fprintf(' -------------- neighbors in mesh --------------------- \n')
neigh_mesh

fprintf(' -------------- neighbors in graph --------------------- \n')
neigh_graph

fprintf(' -------------- type of neighbors ---------------------- \n')
neigh_type
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(a)
colorbar('off')
title('Mesh-grid')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

figure;

subplot(121)
fancy_imagesc(a_index)
colormap(rainbow2_cb(1))
title('Mesh-grid indexes')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
% this is just for visualizing
a_graph = nan(nz,nx);
a_graph(graph2mesh) = (1:n_g2m);
% pcolor does not plot the edges for some weird reason I do not understand
a_graph_plot = [a_graph;nan(1,nx)];
a_graph_plot = [a_graph_plot , nan(nz+1,1)];

% figure;
subplot(122)
fancy_pcolor(a_graph_plot)
colormap(rainbow2_cb(1))
title('Mesh-graph indexes')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(L);
colormap(rainbow2_cb(1))
title('L graph-operator')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
figure;
fancy_pcolor(u_plot)      
colormap(rainbow2_cb(1))
colorbar('off')
title('Mesh-graph field')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
