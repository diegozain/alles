function [ix,iz] = get_izx(izx,nx,nz)
% diego domenzain
% dec 2021
% ------------------------------------------------------------------------------
% get z coordinate
iz = mod(izx,nz);
if (iz==0)
  iz=nz;
end
% ixz = (ix-1)*nz + iz
% get x coordinate
ix = ((izx-iz)/nz) + 1;
end
