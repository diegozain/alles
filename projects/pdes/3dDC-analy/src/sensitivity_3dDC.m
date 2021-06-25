function psi_ = sensitivity_3dDC(X,Z,dx,dz,source_p,source_n,rec_p,rec_n)
% ------------------------------------------------------------------------------
u_p = onelayered_3dDC(1,X,Z,source_p,1);
u_n = onelayered_3dDC(1,X,Z,source_n,1);
u  = u_p - u_n;

u_p = onelayered_3dDC(1,X,Z,rec_p,1);
u_n = onelayered_3dDC(1,X,Z,rec_n,1);
v  = u_p - u_n;
uz=differentiate_plane(u,dz);
ux=differentiate_plane(u.',dx);
ux=ux.';
vz=differentiate_plane(v,dz);
vx=differentiate_plane(v.',dx);
vx=vx.';
psi_ = (ux .* vx) + (uz .* vz);

end