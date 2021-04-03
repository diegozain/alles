function n_g2m = n_g2m_(a,nx,nz)
% diego domenzain
% April 2021 @ Colorado School of Mines
% 
% it is assumed 'a' is a matrix with entries of 0 and 1,
% 0 : a point of no interest
% 1 : a point of interest
% ------------------------------------------------------------------------------
% get the number of nodes in the graph of the mesh 'a'
% ------------------------------------------------------------------------------
b=0;
n_g2m = 0;
for ia = 1:nx*nz
    b=a(ia);
    if b==1
        n_g2m = n_g2m + 1;
    end
end
end