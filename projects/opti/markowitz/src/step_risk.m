function step_ = step_risk(g,k,w,S,E_)

nk=numel(k);
E__ = zeros(nk,1);
for ik=1:nk
  w_= w - k(ik)*g;
  E__(ik) = E_risk(w,S);
end
E__ = [E_;E__];
p = polyfit([0; k],E__,2);
step_ = -p(2)/(2*p(1));
end