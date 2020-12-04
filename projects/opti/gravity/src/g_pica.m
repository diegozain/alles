function [step_x,step_z] = g_pica(d_x,d_z,e_x,e_z,M,Lx,Lz,rho,drho,k)
% diego domenzain
% fall 2017
% Boise State University
% perturb density
rho = rho + drho;
% fwd of perturbed
[ux,uz] = g_fwd(Lx,Lz,rho);
% data of perturbed
d_x_ = M*ux;
d_z_ = M*uz;
% finite difference
d_x_ = d_x_ - d_x;
d_z_ = d_z_ - d_z;
% step size
step_x = k * (d_x_.'*e_x) / (d_x_.'*d_x_);
step_z = k * (d_z_.'*e_z) / (d_z_.'*d_z_);
end