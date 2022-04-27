function n_g2m_3d_(a,nx,ny,nz)
# diego domenzain
# matlab original: August 2021
#
# it is assumed 'a' is a cube matrix with entries of 0 and 1,
# 0 : a point of no interest
# 1 : a point of interest
# ------------------------------------------------------------------------------
# get the number of nodes in the graph of the mesh 'a'
# ------------------------------------------------------------------------------
b=0;
n_g2m = 0;
for ia = 1:nx*ny*nz
    b=a[ia];
    if b==1
        n_g2m = n_g2m + 1;
    end
end
return n_g2m
end
