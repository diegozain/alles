function DX(n,m)
  # ----------------------------------------------------------------------------
  #   x derivative
  # 
  #  on a transposed matrix F of size m by n where m is x and n is z
  # 
  #  Dx = DX(n,m);
  #  dx_F = Dx * F[:]
  # ----------------------------------------------------------------------------
  # second degree accurate,
  # first derivative matrix for
  # 1D discretization on
  # n pts.
  #
  # fwd:             −1.5 	 2 	−0.5
  # ctd:        −0.5 	 0 	 0.5
  # bwd:   0.5 	−2 	 1.5
  # ----------------------------------------------------------------------------
  # fwd
  dx_fwd = [-1.5 2 -0.5];  
  dx_fwd = [dx_fwd zeros(1,n-3)];
  dx_fwd = sparse(dx_fwd);
  # ----------------------------------------------------------------------------
  # centered
  dx_ctd = [-0.5*ones(1,n-2) zeros(1,n-2) 0.5*ones(1,n-2)];
  dx_ctd = dropdims(dx_ctd;dims=1);
  I = 1:n-2;
  I = repeat(transpose(I),1,3);
  I = dropdims(I;dims=1);
  J = [(1:n-2) ; (1+1:n-2+1) ; (1+1+1:n-2+1+1)];
  #
  dx_ctd = sparse(I,J,dx_ctd);
  # ----------------------------------------------------------------------------
  # backward
  dx_bwd = [0.5 -2 1.5];
  dx_bwd = [zeros(1,n-3) dx_bwd];
  dx_bwd = sparse(dx_bwd);
  # ----------------------------------------------------------------------------
  # together
  Dx = [dx_fwd; dx_ctd; dx_bwd];
  # ----------------------------------------------------------------------------
  # make many tiny Dx
  id = sparse(1.0*LinearAlgebra.I, m, m); #speye(m);
  Dx = kron(id,Dx);
  return Dx
end