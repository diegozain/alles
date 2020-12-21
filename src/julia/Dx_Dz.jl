function Dx_Dz(n,m)
  # diego domenzain
  # fall 2018 @ BSU
  # ----------------------------------------------------------------------------
  # first degree derivative, 
  # second degree accurate.
  # 
  # fwd:             −1.5 	 2 	−0.5
  # ctd:        −0.5 	 0 	 0.5
  # bwd:   0.5 	−2 	 1.5
  #
  # Dx: vertical derivative of a matrix
  # Dz: horizontal derivative of a matrix.
  # they are for transposed matricies. dont whine.
  # ----------------------------------------------------------------------------
  # example:
  # [nz,nx] = size(a);
  # [Dz,Dx] = Dx_Dz(nz,nx);
  # ----------------------------------------------------------------------------
  Dx = DX(n,m);
  Dz = DZ(n,m);
  return Dx,Dz
end