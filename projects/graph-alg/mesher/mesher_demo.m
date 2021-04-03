clear
close all
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
% ------------------------------------------------------------------------------
% since we are assuming a 2D mesh-grid, max # of neighbors in the mesh is 4
n_neigh_max = 4;
% ------------------------------------------------------------------------------
% setup a simple example
a = [0 0 0 1 1 0 0; 0 1 1 1 1 0 0; 1 1 1 1 1 1 1; 1 1 0 0 0 1 1];
[nz,nx]=size(a);
% ------------------------------------------------------------------------------
% this is only for easy reference:
a_index = 1:(nx*nz);
a_index = reshape(a_index,[nz,nx]);
% ------------------------------------------------------------------------------
% get the number of nodes in the graph 
b=0;
n_g2m = 0;
for ia = 1:nx*nz
    b=a(ia);
    if b==1
        n_g2m = n_g2m + 1;
    end
end
% ------------------------------------------------------------------------------
% make two dictionaries,
% graph2mesh : indexes are graph nodes, entries are mesh nodes
% mesh2graph : indexes are mesh nodes, entries are graph nodes
graph2mesh = zeros(n_g2m,1);
mesh2graph = zeros(nz*nx,1);
i_g2m = 0;
for ia = 1:nx*nz
    b=a(ia);
    if b==1
        i_g2m = i_g2m + 1;
        graph2mesh(i_g2m) = ia;
        mesh2graph(ia) = i_g2m;
    end
end
% ------------------------------------------------------------------------------
% each node has neighbors.
% neigh_mesh : row indexes are graph nodes.
%              row entries are neighbors of that node, in the mesh. 
neigh_mesh = zeros(n_g2m,n_neigh_max);

i_up = 0;
i_do = 0;
i_ri = 0;
i_le = 0;

for i_g2m = 1:n_g2m
    
    i_ri = graph2mesh(i_g2m) + nz;
    i_up = graph2mesh(i_g2m) - 1;
    i_le = graph2mesh(i_g2m) - nz;
    i_do = graph2mesh(i_g2m) + 1;
    
    % left edge
    if (graph2mesh(i_g2m)<=nz)
        if (a(i_ri)==1)
            neigh_mesh(i_g2m,1) = i_ri;
        end
        
        if i_up>=1
            if (a(i_up)==1)
                neigh_mesh(i_g2m,2) = i_up;
            end
        end
        if i_do<=nz
            if (a(i_do)==1)
                neigh_mesh(i_g2m,4) = i_do;
            end
        end
    % right edge
    elseif (graph2mesh(i_g2m)>nz*(nx-1))
        if (a(i_le)==1)
            neigh_mesh(i_g2m,3) = i_le;
        end
        
        if i_up>=nz*(nx-1)+1
            if (a(i_up)==1)
                neigh_mesh(i_g2m,2) = i_up;
            end
        end
        if i_do<=(nz*nx)
            if (a(i_do)==1)
                neigh_mesh(i_g2m,4) = i_do;
            end
        end
    % bottom edge
    elseif (mod(graph2mesh(i_g2m),nz)==0)
        if (a(i_up)==1)
            neigh_mesh(i_g2m,2) = i_up;
        end
        
        if i_ri<=(nz*nx)
            if (a(i_ri)==1)
                neigh_mesh(i_g2m,1) = i_ri;
            end
        end
        if i_le>=nz
            if (a(i_le)==1)
                neigh_mesh(i_g2m,3) = i_le;
            end
        end
    % top edge
    elseif (mod(graph2mesh(i_g2m),nz)==1)
        if (a(i_do)==1)
            neigh_mesh(i_g2m,4) = i_do;
        end
        
        if i_ri<=(nz*(nx-1)+1)
            if (a(i_ri)==1)
                neigh_mesh(i_g2m,1) = i_ri;
            end
        end
        if i_le>=1
            if (a(i_le)==1)
                neigh_mesh(i_g2m,3) = i_le;
            end
        end
    % inner nodes
    else
        if (a(i_ri)==1)
            neigh_mesh(i_g2m,1) = i_ri;
        end
        if (a(i_up)==1)
            neigh_mesh(i_g2m,2) = i_up;
        end
        if (a(i_le)==1)
            neigh_mesh(i_g2m,3) = i_le;
        end
        if (a(i_do)==1)
            neigh_mesh(i_g2m,4) = i_do;
        end
    end
end
% ------------------------------------------------------------------------------
% each node has neighbors.
% neigh_graph : row indexes are graph nodes.
%               row entries are neighbors of that node, in the graph.
neigh_graph = zeros(n_g2m,n_neigh_max);
for i_g2m = 1:n_g2m
    for i_nei = 1:n_neigh_max
        if (neigh_mesh(i_g2m,i_nei) ~= 0)
            neigh_graph(i_g2m,i_nei) = mesh2graph(neigh_mesh(i_g2m,i_nei));
        end
    end
end
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
neigh_type = zeros(n_g2m,n_neigh_max);

inner =  1;
neuma = -1;

i_up = 0;
i_do = 0;
i_ri = 0;
i_le = 0;

for i_g2m = 1:n_g2m
    
    i_ri = graph2mesh(i_g2m) + nz;
    i_up = graph2mesh(i_g2m) - 1;
    i_le = graph2mesh(i_g2m) - nz;
    i_do = graph2mesh(i_g2m) + 1;
    
    % left edge
    if (graph2mesh(i_g2m)<=nz)
        if (a(i_ri)==1)
            neigh_type(i_g2m,1) = inner;
        else
            neigh_type(i_g2m,1) = neuma;
        end
        
        if i_up>=1
            if (a(i_up)==1)
                neigh_type(i_g2m,2) = inner;
            else
                neigh_type(i_g2m,2) = neuma;
            end
        end
        if i_do<=nz
            if (a(i_do)==1)
                neigh_type(i_g2m,4) = inner;
            else
                neigh_type(i_g2m,4) = neuma;
            end
        end
    % right edge
    elseif (graph2mesh(i_g2m)>nz*(nx-1))
        if (a(i_le)==1)
            neigh_type(i_g2m,3) = inner;
        else
            neigh_type(i_g2m,3) = neuma;
        end
        
        if i_up>=nz*(nx-1)+1
            if (a(i_up)==1)
                neigh_type(i_g2m,2) = inner;
            else
                neigh_type(i_g2m,2) = neuma;
            end
        end
        if i_do<=(nz*nx)
            if (a(i_do)==1)
                neigh_type(i_g2m,4) = inner;
            else
                neigh_type(i_g2m,4) = neuma;
            end
        end
    % bottom edge
    elseif (mod(graph2mesh(i_g2m),nz)==0)
        if (a(i_up)==1)
            neigh_type(i_g2m,2) = inner;
        else
            neigh_type(i_g2m,2) = neuma;
        end
        
        if i_ri<=(nz*nx)
            if (a(i_ri)==1)
                neigh_type(i_g2m,1) = inner;
            else
                neigh_type(i_g2m,1) = neuma;
            end
        end
        if i_le>=nz
            if (a(i_le)==1)
                neigh_type(i_g2m,3) = inner;
            else
                neigh_type(i_g2m,3) = neuma;
            end
        end
    % top edge
    elseif (mod(graph2mesh(i_g2m),nz)==1)
        if (a(i_do)==1)
            neigh_type(i_g2m,4) = inner;
        else
            neigh_type(i_g2m,4) = neuma;
        end
        
        if i_ri<=(nz*(nx-1)+1)
            if (a(i_ri)==1)
                neigh_type(i_g2m,1) = inner;
            else
                neigh_type(i_g2m,1) = neuma;
            end
        end
        if i_le>=1
            if (a(i_le)==1)
                neigh_type(i_g2m,3) = inner;
            else
                neigh_type(i_g2m,3) = neuma;
            end
        end
    % inner nodes
    else
        if (a(i_ri)==1)
            neigh_type(i_g2m,1) = inner;
        else
            neigh_type(i_g2m,1) = neuma;
        end
        if (a(i_up)==1)
            neigh_type(i_g2m,2) = inner;
        else
            neigh_type(i_g2m,2) = neuma;
        end
        if (a(i_le)==1)
            neigh_type(i_g2m,3) = inner;
        else
            neigh_type(i_g2m,3) = neuma;
        end
        if (a(i_do)==1)
            neigh_type(i_g2m,4) = inner;
        else
            neigh_type(i_g2m,4) = neuma;
        end
    end
end
% ------------------------------------------------------------------------------
fprintf(' mesh under consideration\n')
a

fprintf(' indexes of previous mesh\n')
a_index

fprintf(' graph -> mesh\n')
graph2mesh

fprintf(' mesh -> graph\n')
mesh2graph

fprintf(' neighbors in mesh\n')
neigh_mesh

fprintf(' neighbors in graph\n')
neigh_graph

fprintf(' type of neighbors\n')
neigh_type
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(a)
colorbar('off')
title('Mesh-grid')
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
