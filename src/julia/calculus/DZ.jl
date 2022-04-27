function DZ(n,m)
  # ----------------------------------------------------------------------------
  #   z derivative
  # 
  #  on a transposed matrix F of size m by n where m is x and n is z
  # 
  #  Dz = DZ(n,m);
  #  dz_F = Dz * F[:]
  # ----------------------------------------------------------------------------
  # first degree derivative, 
  # second degree accurate.
  # 
  # fwd:             −1.5 	 2 	−0.5
  # ctd:        −0.5 	 0 	 0.5
  # bwd:   0.5 	−2 	 1.5
  # ----------------------------------------------------------------------------
  # fwd
  dz_fwd = [-1.5*ones(1,n) 2*ones(1,n) -0.5*ones(1,n) zeros(1,n*m-3*n)];
  dz_fwd = dropdims(dz_fwd;dims=1);
  I = 1:n;
  I = repeat(transpose(I),1,m);
  I = dropdims(I;dims=1);
  J = 1:n*m;
  dz_fwd = sparse(I,J,dz_fwd);
  # ----------------------------------------------------------------------------
  # center
  dz_ctd = [-0.5*ones(1,n*m-2*n) zeros(1,n*m-2*n) 0.5*ones(1,n*m-2*n)];
  dz_ctd = dropdims(dz_ctd;dims=1);
  I = 1:n*m-2*n;
  I = repeat(transpose(I),1,3);
  I = dropdims(I;dims=1);
  J = [(1:n*m-2*n) ; (n+1:n*m-2*n+n) ; (n+n+1:n*m-2*n+n+n)];
  dz_ctd = sparse(I,J,dz_ctd);
  # ----------------------------------------------------------------------------
  # bwd
  dz_bwd = [zeros(1,n*m-3*n) 0.5*ones(1,n) -2*ones(1,n) 1.5*ones(1,n)];
  dz_bwd = dropdims(dz_bwd;dims=1);
  I = 1:n;
  I = repeat(transpose(I),1,m);
  I = dropdims(I;dims=1);
  J = 1:n*m;
  dz_bwd = sparse(I,J,dz_bwd);
  # ----------------------------------------------------------------------------
  # together
  Dz = [dz_fwd; dz_ctd; dz_bwd];
  #
  return Dz
end