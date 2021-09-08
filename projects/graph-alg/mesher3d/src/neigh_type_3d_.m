function neigh_type = neigh_type_3d_(a,nx,ny,nz,n_g2m,graph2mesh)
% diego domenzain
% August 2021
%
% ------------------------------------------------------------------------------
% each node has a special type in a mesh-grid.
% neighbors that are (in the mesh):
% zero, non-zero, and next to the limits of the mesh.
%
% In the case the mesh is a chunk of the earth, and we are doing this for a PDE,
% this translates to:
%
% neighbors that are non-zero = the pde (inner-node)
% neighbors that are zero = neumann bc
% neighbors that are next to the limits of the mesh = robin bc
%
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
%
% we define : (type,BC) = (1,inner) (-1,neumann) (0,robin)
% ------------------------------------------------------------------------------
neigh_type = zeros(n_g2m,6,'int32');

inner =  1;
neuma = -1;

iyxz = 0;
inei = 0;

i_up = 0;
i_do = 0;
i_ri = 0;
i_le = 0;
i_fw = 0;
i_bw = 0;

for i_g2m = 1:n_g2m
    iyxz = graph2mesh(i_g2m);
    % get array index coordinates of all 6 neighbors
    i_ri = iyxz + nz;
    i_up = iyxz - 1;
    i_le = iyxz - nz;
    i_do = iyxz + 1;
    i_fr = iyxz - nz*nx;
    i_ba = iyxz + nz*nx;
    % --------------------------------------------------------------------------
    % get cube index coordinates of all 6 neighbors
    [ix,iy,iz] = get_ixyz(iyxz,nx,ny,nz);
    % --------------------------------------------------------------------------
    % now access matrix cube 'a'
    % --------------------------------------------------------------------------
    % right neighbor
    inei = ix + 1;
    % access possible domain only
    if (inei <= nx)
      if (a(i_ri) == 1)
          neigh_type(i_g2m,1) = inner;
      else
          neigh_type(i_g2m,1) = neuma;
      end
    end
    % up neighbor
    inei = iz - 1;
    % access possible domain only
    if (inei >= 1)
      if (a(i_up) == 1)
          neigh_type(i_g2m,2) = inner;
      else
          neigh_type(i_g2m,2) = neuma;
      end
    end
    % left neighbor
    inei = ix - 1;
    % access possible domain only
    if (inei >= 1)
      if (a(i_le) == 1)
          neigh_type(i_g2m,3) = inner;
      else
          neigh_type(i_g2m,3) = neuma;
      end
    end
    % down neighbor
    inei = iz + 1;
    % access possible domain only
    if (inei <= nz)
      if (a(i_do) == 1)
          neigh_type(i_g2m,4) = inner;
      else
          neigh_type(i_g2m,4) = neuma;
      end
    end
    % front neighbor
    inei = iy - 1;
    % access possible domain only
    if (inei >= 1)
      if (a(i_fr) == 1)
          neigh_type(i_g2m,5) = inner;
      else
          neigh_type(i_g2m,5) = neuma;
      end
    end
    % back neighbor
    inei = iy + 1;
    % access possible domain only
    if (inei <= ny)
      if (a(i_ba) == 1)
          neigh_type(i_g2m,6) = inner;
      else
          neigh_type(i_g2m,6) = neuma;
      end
    end
end
end
