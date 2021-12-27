function neigh_mesh = neigh_mesh_(a,nx,nz,n_g2m,graph2mesh)
% diego domenzain
% April 2021 @ Colorado School of Mines
%
% it is assumed 'a' is a matrix with entries of 0 and 1,
% 0 : a point of no interest
% 1 : a point of interest
% ------------------------------------------------------------------------------
% neigh_mesh : row indexes are graph nodes.
%              row entries are neighbors of that node, in the mesh.
%
%      2
%      |
% 3 -- i -- 1
%      |
%      4
% ------------------------------------------------------------------------------
neigh_mesh = zeros(n_g2m,4,'uint32');

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
end
