function p = p2ij(p,nx)

p_i=mod(p,nx); p_i(p_i==0)=nx;
p_j=((p-p_i)/nx)+1;
p = [full(p_i) full(p_j)];

end