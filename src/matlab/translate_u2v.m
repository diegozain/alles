function f = translate_u2v(u,v)
% diego domenzain
% oct 2022
% ------------------------------------------------------------------------------
% find translator function f that maps,
%
%                f:u ⟶ v
%                v = f(u)
%
% u & v must be arrays of indexes (ℕ)
% f will also be an array of indexes (ℕ)
% ------------------------------------------------------------------------------
% example
%
% u =[1; 2; 4; 7; 8; 10; 13; 14; 16];
% v =[1; 2; 3; 4; 5;  6;  7;  8;  9];
%
% then,
%          f = [1 ; 2; 0; 3; 0; 0; 4; 5; 0; 6; 0; 0; 7; 8; 0; 9];
%
% ------------------------------------------------------------------------------
nu = numel(u);
nf = max(u);
f  = zeros(nu,1,'uint32');
for iu=1:nu
 f(u(iu))=v(iu);
end
end
