function Lw = L_w(p,s_w)

har = @(u,v) (2*u*v) / (u+v);

a = har( p(1) , p(1) );
b = har( p(1) , p(2) );
c = har( p(2) , p(2) );

Lw = [a-b b;...
      -b c-b];

end
