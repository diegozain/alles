function integrate_cube_(v,dt)
% diego domenzain
% spring 2021 @ CSM
% ------------------------------------------------------------------------------
% 
% given a time cube 'v' sampled at dt, who is its antiderivative?
% 
% clearly, integral(v,dt) is unique up to a constant. 
% here we assume the constant of integration is zero.
% 
% tight memory requirements: no extra copy of 'v' is stored.
% ------------------------------------------------------------------------------
% this method follows the integration scheme in integrate_line.m,
% and does so as a void method. 
% 
% the idea of this method is explained in detail in Idt.m
% ------------------------------------------------------------------------------
a = 0.5*dt;
b = 2*dt;

[nz,nx,nt] = size(v.u);

v_ = zeros(nz,nx);
v__= zeros(nz,nx);
% ------------------------------------------------------------------------------
v_= v.u(:,:,3);
v.u(:,:,3) = b * v.u(:,:,2);
v.u(:,:,2) = a * (v.u(:,:,1) + v.u(:,:,2));
for it=4:nt
 v__  = v.u(:,:,it);
 v.u(:,:,it)= v.u(:,:,it-2) + b*v_;
 v_   = v__;
end
v.u(:,:,1) = zeros(nz,nx);
end