function integrate_line_(v,dt)
% diego domenzain
% spring 2021 @ CSM
% ------------------------------------------------------------------------------
% 
% given a time series 'y' sampled at dt, who is its antiderivative?
% 
% clearly, integral(y,dt) is unique up to a constant. 
% we assume the constant of integration is zero.
% ------------------------------------------------------------------------------
% this code does the idea of Idt.m but no matrix storage is needed.
% if you ** dont ** want to visualize the matrix Idt_, use this method.
% ------------------------------------------------------------------------------
a = 0.5*dt;
b = 2*dt;

nt= numel(v.u);

y_ =0;
y__=0;
% ------------------------------------------------------------------------------
y_ = v.u(3);
v.u(3) = b * v.u(2);
v.u(2) = a * (v.u(1) + v.u(2));
for it=4:nt
 y__  = v.u(it);
 v.u(it)= v.u(it-2) + b*y_;
 y_   = y__;
end
v.u(1) = 0;
end
