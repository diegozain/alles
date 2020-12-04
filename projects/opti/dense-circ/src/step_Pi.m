function step_ = step_Pi(r,g,lam,expo,Pi_,k,ir)

nk=numel(k);
Pi__ = zeros(nk,1);
for ik=1:nk
  r_=r;
  r_(ir,:)= r(ir,:) - k(ik)*g;
  Pi__(ik)  = Pi(r_,lam,expo,ir);
end
Pi__ = [Pi_;Pi__];
p = polyfit([0; k],Pi__,2);
step_ = -p(2)/(2*p(1));
end