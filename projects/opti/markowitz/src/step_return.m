function step_ = step_return(g,k,w,r,rho_,E_)

nk=numel(k);
E__ = zeros(nk,1);
for ik=1:nk
  w_= w - k(ik)*g;
  E__(ik) = E_return(w,r,rho_);
end
E__ = [E_;E__];
p = polyfit([0; k],E__,2);
step_ = -p(2)/(2*p(1));
end