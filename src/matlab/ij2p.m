function p = ij2p(p,nx)
ip = p(1);
jp = p(2);
p = ip + (jp-1)*nx;
end