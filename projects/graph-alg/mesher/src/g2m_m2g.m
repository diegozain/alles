function [graph2mesh,mesh2graph] = g2m_m2g(a,nx,nz,n_g2m)
% diego domenzain
% April 2021 @ Colorado School of Mines
% 
% it is assumed 'a' is a matrix with entries of 0 and 1,
% 0 : a point of no interest
% 1 : a point of interest
% ------------------------------------------------------------------------------
% make two dictionaries,
% graph2mesh : indexes are graph nodes, entries are mesh nodes
% mesh2graph : indexes are mesh nodes, entries are graph nodes
% ------------------------------------------------------------------------------
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
end