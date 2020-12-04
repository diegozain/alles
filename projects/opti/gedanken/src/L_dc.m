function Ldc = L_dc(p,s_dc)

har = @(u,v) (2*u*v) / (u+v);

a = har( p(1) , p(1) );
b = har( p(1) , p(2) );
c = har( p(2) , p(2) );

Ldc = [a-b c;...
      -b c-b];

end
