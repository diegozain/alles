function [g,S] = Jte_dc(u,s,p,e,M,L)

sa = M' * e;
a = adj(L,sa);
S = S_dc(u,p,s);
g = S*a;

end
