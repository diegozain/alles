function [ux,uz] = g_fwd(Lx,Lz,rho)
% diego domenzain
% fall 2017
% Boise State University
uz = Lz*rho;
ux = Lx*rho;
end