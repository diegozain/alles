function w = Jg_dc(u,s,p,g,M,L)

S = S_dc(u,p,s);
w = S' * g;
w = L \ w;
w = M * w;

end
