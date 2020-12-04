function step_ = step_Mi(r,g,m,Mi_,k,ir)

nk=numel(k);
Mi__ = zeros(nk,1);
for ik=1:nk
  r_=r;
  r_(ir,:)= r(ir,:) - k(ik)*g;
  Mi__(ik)  = Mi(r_,m,ir);
end
Mi__ = [Mi_;Mi__];
p = polyfit([0; k],Mi__,2);
step_ = -p(2)/(2*p(1));
end