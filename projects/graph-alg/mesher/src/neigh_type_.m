function neigh_type = neigh_type_(a,nx,nz,n_g2m,graph2mesh)
% diego domenzain
% April 2021 @ Colorado School of Mines
% 
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
neigh_type = zeros(n_g2m,4,'int32');

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
end