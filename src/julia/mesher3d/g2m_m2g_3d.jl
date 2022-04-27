function g2m_m2g_3d(a,nx,ny,nz,n_g2m)
# diego domenzain
# matlab original: August 2021
#
# it is assumed 'a' is a cube matrix with entries of 0 and 1,
# 0 : a point of no interest
# 1 : a point of interest
# ------------------------------------------------------------------------------
# make two dictionaries,
# graph2mesh : indexes are graph nodes, entries are mesh nodes
# mesh2graph : indexes are mesh nodes, entries are graph nodes
# ------------------------------------------------------------------------------
graph2mesh = zeros(Int32,n_g2m,1);
mesh2graph = zeros(Int32,nx*ny*nz,1);
i_g2m = 0;
for ia = 1:nx*ny*nz
    b=a[ia];
    if b==1
        i_g2m = i_g2m + 1;
        graph2mesh[i_g2m] = ia;
        mesh2graph[ia] = i_g2m;
    end
end
return graph2mesh, mesh2graph
end
