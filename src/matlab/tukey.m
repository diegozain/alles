function tuk = tukey(nt,r)

x=linspace(0,1,nt).';
tuk=zeros(nt,1);
for it=1:nt
 if (0 <= x(it)) &&  (x(it) < r*0.5)
  tuk(it,1) = 0.5*(1+cos((2*pi/r) * (x(it) - r*0.5) ));
 end
 if (r*0.5 <= x(it)) && (x(it) < (1-r*0.5))
  tuk(it,1) = 1;
 end
 if ((1-r*0.5) <= x(it)) && (x(it) <= 1)
  tuk(it,1) = 0.5*(1 + cos( (2*pi/r) * ( x(it) - 1 + r*0.5 ) ));
 end
end
end