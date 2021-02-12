function [Lx,Lz] = g_L(x,z)
% diego domenzain
% fall 2017
% Boise State University
% -- in
% spacial discretizations.
% -- out
% linear operator of 'sum of forces'.
% two integrals, one for x (Lx) one for z (Lz).
% --
% boundary conditions are Dirichlet, 
% only good for examples, not for field data.
dx = x(2)-x(1) ;
dz = z(2)-z(1) ;

nx = numel(x);
nz = numel(z);

nu = nx*nz;
Lz = zeros(nu,nu);
Lx = zeros(nu,nu);

X = repmat(x,[1,nz]);
Z = repmat(z,[1,nx]).';

for ix = 1:nx;
	for iz = 1:nz;
		Xi = X-x( ix );
		Zi = Z-z( iz );
		R = sqrt(Xi.^2 + Zi.^2);
		% this accounts for R(ix,iz) = 0 (dont want no NaN below)
		R( ix,iz ) = 1;
		
		% abs because there are no gravity dipoles
		Lzi = abs(Zi)./(R.^3);
		Lxi = abs(Xi)./(R.^3);

		p = ij2p([ix,iz],nx);
		Lz(p,:) = (dx*dz)*Lzi(:).';
		Lx(p,:) = (dx*dz)*Lxi(:).';

		% this accounts for R(ix,iz) = 0, we still want zero there (but no NaN)
		Lz(p,p)=0;
		Lx(p,p)=0;
	end
end
end