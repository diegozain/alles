function [u_w,e_w,Lw] = fwd_w(p,s_w,Mw,d_w_o) 

Lw = L_w(p,s_w);

u_w = Lw \ s_w;
d = Mw * u_w;
e_w = d - d_w_o;

end
