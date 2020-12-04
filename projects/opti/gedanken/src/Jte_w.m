function [g,S] = Jte_w(u,s,p,e,M,L)

sa = M' * e;
a = adj(L,sa);
S = S_w(u,p,s);
g = S*a;

end
