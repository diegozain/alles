function v = haar1d_inv(v)
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
s1_2 = sqrt(1.0/2);
n = numel(v);
h = sqrt(n);

v = v*h;

n_=log2(n);

for d_=(n_-1):-1:0
 d=2^d_;
 for i_=0:(2*d):(n-1)
  x = 0.0; y = 0.0;
  x = s1_2 * (v(i_+1) + v(i_+d+1) );
  y = s1_2 * (v(i_+1) - v(i_+d+1) );
  
  v(i_+1)  = x;
  v(i_+d+1)= y;
 end
end
end
% ------------------------------------------------------------------------------
% % TEST
% % b_ should be equal to b
% % a_ should be equal to a
% 
% % coefficients of "block" basis
% a = [1.0000 ; 0.5000 ; 0.3333 ; 0.2500 ; 0.2000 ; 0.1667 ; 0.1429 ; 0.1250];
% % coefficients of "weird" basis
% b = [0.3397 ; 0.1250 ; 0.1620 ; 0.0208 ; 0.1811 ; 0.0083 ; 0.0175 ; 0.0045];
% 
% b_= haar1d_fwd(a);
% a_= haar1d_inv(b);
% ------------------------------------------------------------------------------

