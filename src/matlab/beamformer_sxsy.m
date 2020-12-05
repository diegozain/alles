function B = beamformer_sxsy(fo,r,d_,f,sX,sY)

nfo = numel(fo);
[nsy,nsx] = size(sX);
B = zeros( nsy, nsx , nfo);

for i=1:nfo
B(:,:,i) = beamformer(fo(i),r,d_,f,sX,sY);
end
B = sum(B,3);
end