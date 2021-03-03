function v = haar2d_inv(v)
% ------------------------------------------------------------------------------
% diego domenzain 
% @ CSM spring 2021
%
% matlab implementation of the inverse 1d haar basis.
% see src/c/haar.c for the faster better C implementation.
%
% from the book
% Programming Projects in C for Students of Engineering, Science, and Mathematics
% by Rouben Rostamian.
%
% only takes v of size a power of 2.
% ------------------------------------------------------------------------------
% since we're in matlab, we can abuse it a bit and not worry about column/row 
[nr,nc]=size(v);
ic=1:nc;
v(:,ic) = haar1d_inv(v(:,ic));
ir=1:nr;
v(ir,:) = haar1d_inv(v(ir,:));
end
% ------------------------------------------------------------------------------
% % TEST
% ??
% ------------------------------------------------------------------------------
