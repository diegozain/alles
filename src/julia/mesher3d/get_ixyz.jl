function get_ixyz(iyxz,nx,ny,nz)
# diego domenzain
# matlab original: August 2021
# ------------------------------------------------------------------------------
# get z coordinate
iz = mod(iyxz,nz);
if (iz==0)
  iz=nz;
end
# iyxz = (iy-1)*nx*nz + (ix-1)*nz + iz  ... (*)
ixz = mod(iyxz,nx*nz);
if (ixz==0)
  ixz=nx*nz;
end
# ixz = (ix-1)*nz + iz from (*)
# get x coordinate
ix = ((ixz-iz)/nz) + 1;
# get y coordinate from (*)
iy = ((iyxz-ixz)/(nx*nz)) + 1;
#
ixyz_ = [ix,iy,iz];
ixyz_=convert(Array{Int32},ixyz_)
return ixyz_
end
