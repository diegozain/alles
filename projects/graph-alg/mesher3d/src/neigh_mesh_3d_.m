function neigh_mesh = neigh_mesh_3d_(a,nx,ny,nz,n_g2m,graph2mesh)
% diego domenzain
% August 2021
%
% it is assumed 'a' is a cube matrix with entries of 0 and 1,
% 0 : a point of no interest
% 1 : a point of interest
% ------------------------------------------------------------------------------
% neigh_mesh : row indexes are graph nodes.
%              row entries are neighbors of that node, in the mesh.
% ------------------------------------------------------------------------------
neigh_mesh = zeros(n_g2m,6,'uint32');

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
    % get z coordinate
    iz = mod(iyxz,nz);
    if (iz==0)
      iz=nz;
    end
    % iyxz = (iy-1)*nx*nz + (ix-1)*nz + iz  ... (*)
    ixz = mod(iyxz,nx*nz);
    if (ixz==0)
      ixz=nx*nz;
    end
    % ixz = (ix-1)*nz + iz from (*)
    % get x coordinate
    ix = ((ixz-iz)/nz) + 1;
    % get y coordinate from (*)
    iy = ((iyxz-ixz)/(nx*nz)) + 1;
    % --------------------------------------------------------------------------
    % now access matrix cube 'a'
    % --------------------------------------------------------------------------
    % right neighbor
    inei = ix + 1;
    % access possible domain only
    if (inei <= nx)
      if (a(i_ri) == 1)
          neigh_mesh(i_g2m,1) = i_ri;
      end
    end
    % up neighbor
    inei = iz - 1;
    % access possible domain only
    if (inei >= 1)
      if (a(i_up) == 1)
          neigh_mesh(i_g2m,2) = i_up;
      end
    end
    % left neighbor
    inei = ix - 1;
    % access possible domain only
    if (inei >= 1)
      if (a(i_le) == 1)
          neigh_mesh(i_g2m,3) = i_le;
      end
    end
    % down neighbor
    inei = iz + 1;
    % access possible domain only
    if (inei <= nz)
      if (a(i_do) == 1)
          neigh_mesh(i_g2m,4) = i_do;
      end
    end
    % front neighbor
    inei = iy - 1;
    % access possible domain only
    if (inei >= 1)
      if (a(i_fr) == 1)
          neigh_mesh(i_g2m,5) = i_fr;
      end
    end
    % back neighbor
    inei = iy + 1;
    % access possible domain only
    if (inei <= ny)
      if (a(i_ba) == 1)
          neigh_mesh(i_g2m,6) = i_ba;
      end
    end
end
end
