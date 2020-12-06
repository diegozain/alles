function b = beamformer(fo,r,d_,f,sX,sY)
% b -> COMPLEX beamform
%
% fo <- desired freq
% r  <- receiver coordinates ( nr x 2 )
% d_  <- data ( nf x nr ) in frequency domain
% f <- frequency vector ( fftshifted & centered )
% [sX,sY] <- meshgrid of x and y slowness.
%
%
% b( wo ) = sum_r { d_(r,wo) * exp( i*wo * r*s ) }
%
%
% b will be a matrix of size
% ( nsy x nsx ) where 
% nsx is # of columns in sX
% nsy is # of rows in sY
%
b = zeros( size(sX) );

% find real # fo in discretized vector f.
% returns index,
% (or if fo is a vector, returns indicies).
%
ifo = binning(f,fo);

% now beamform loop on receivers 
% because dylan hasn't taught me
% matrix algorithm for beamforming :(
%
[~,nr] = size( d_ );
nwo = numel( fo );
for jr=1:nr
	% t_shift is a matrix
	% of size ( nsy x nsx )
	%
	t_shift = ( r(jr,1)*sX + r(jr,2)*sY );
	b = b + d_(ifo,jr) * exp( t_shift * 1i*2*pi*fo );
end
b = b/nr;
end