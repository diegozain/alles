function S = S_w(u,p,s)

% A' = d_p1(L) * field
% B' = d_p2(L) * field
%
% A and B are the matrices of 
% grad_p(u) = [A , B] 
% when u = Ls.
%
% ... although this matrix S is used 
% when Lu = s.

dhar_u = @(u,v) (2*v*v) / (u+v)^2;

b1 = dhar_u( p(1),p(2) );
b2 = dhar_u( p(2),p(1) );

A = [1-b1,  b1; ...
      -b1, -b1];

B = [-b2, b2; ...
      -b2, 1-b2];

A = (A * u)';
B = (B * u)';

S = - [A ; B];

end
