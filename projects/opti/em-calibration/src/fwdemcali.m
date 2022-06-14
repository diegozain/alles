function data = fwdemcali(param,s)
% diego domenzain
% jun 2022
% ------------------------------------------------------------------------------
l = param(1);
c = param(2);
r = param(3);
a = param(4);
b = param(5);
% ------------------------------------------------------------------------------
data = ( l*c*s.^2 + (r*c + l/a)*s + 1 + r/a ) ./ ( l*c*s.^2 + (r*c + l/b)*s + 1 + r/b );
end
