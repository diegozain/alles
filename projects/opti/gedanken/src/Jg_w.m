function w = Jg_w(u,s,p,g,M,L)

S = S_w(u,p,s);
w = S' * g;
w = L \ w;
w = M * w;

end
