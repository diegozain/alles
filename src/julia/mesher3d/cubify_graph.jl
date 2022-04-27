function cubify_graph(phi_,graph2mesh,nx,ny,nz,n_g2m)
# diego domenzain
# matlab original: winter 2021
# ------------------------------------------------------------------------------
# get phi_ in 🍇 and put it in 🎲
# ------------------------------------------------------------------------------
phi3d_=Array{Float64}(undef,ny,nx,nz);
phi3d_[:].=NaN;
for i_g2m=1:n_g2m
  # get x,y,z coordinate
  iyxz = graph2mesh[i_g2m];
  ixyz_= get_ixyz(iyxz,nx,ny,nz);

  phi3d_[ixyz_[2],ixyz_[1],ixyz_[3]]= phi_[i_g2m];
end
return phi3d_
end
