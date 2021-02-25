function y = integrate_(y,dt)
% diego domenzain
% spring 2021 @ CSM
% ------------------------------------------------------------------------------
% 
% given a time series 'y' sampled at dt, who is its antiderivative?
% 
% clearly, integral(y) is unique up to a constant. 
% here we assume the constant of integration is zero.
% ------------------------------------------------------------------------------
% this code does the idea of Idt.m but no matrix storage is needed.
% if you ** dont ** want to visualize the matrix Idt_, use this method.
% ------------------------------------------------------------------------------
a = 0.5*dt;
b = 2*dt;
nt= numel(y);
% ------------------------------------------------------------------------------
y_ = y(3);
y(3) = b * y(2);
y(2) = a * (y(1) + y(2));
for it=4:nt
 y__  = y(it);
 y(it)= y(it-2) + b*y_;
 y_   = y__;
end
y(1) = 0;
end