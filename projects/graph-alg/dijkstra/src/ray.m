function p = ray(s,r,P)
  
p = [];
while P(r) > 0
  p = [r;p];
  r = P(r);
end
p = [s;p];
end