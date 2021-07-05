function psi_ = sensitivity_3dDC(X,Z,dx,dz,source_p,source_n,rec_p,rec_n)
% diego domenzain @ aarhus uni, summer 2021
% 
% compute sensitivities for one abmn pair in a homogeneous background
% a la Spitzer.
% ------------------------------------------------------------------------------
% φ
u_p = onelayered_3dDC(1,X,Z,source_p,1);
u_n = onelayered_3dDC(1,X,Z,source_n,1);
u  = u_p - u_n;
% ------------------------------------------------------------------------------
% v
nrecs = size(rec_p,1);
v = zeros(size(u));
for irecs = 1:nrecs
 u_p = onelayered_3dDC(1,X,Z,rec_p(irecs,:),1);
 u_n = onelayered_3dDC(1,X,Z,rec_n(irecs,:),1);
 v_  = u_p - u_n;
 v = v + v_;
end
% ------------------------------------------------------------------------------
% ∇φ & ∇v
uz=differentiate_plane(u,dz);
ux=differentiate_plane(u.',dx);
ux=ux.';
vz=differentiate_plane(v,dz);
vx=differentiate_plane(v.',dx);
vx=vx.';
% ------------------------------------------------------------------------------
% ∇φ ⋅ ∇v
psi_ = (ux .* vx) + (uz .* vz);
end