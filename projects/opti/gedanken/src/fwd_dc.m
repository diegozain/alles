function [u_dc,e_dc,Ldc] = fwd_dc(p,s_dc,Mdc,d_dc_o) 

Ldc = L_dc(p,s_dc);

u_dc = Ldc \ s_dc;
d = Mdc * u_dc;
e_dc = d - d_dc_o;

end
