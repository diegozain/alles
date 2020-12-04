function step_ = step_Ei(r,g,d,Ei_,k,ir)

nk=numel(k);
Ei__ = zeros(nk,1);
for ik=1:nk
  r_=r;
  r_(ir,:)= r(ir,:) - k(ik)*g;
  Ei__(ik)  = Ei(r_,d,ir);
end
Ei__ = [Ei_;Ei__];
p = polyfit([0; k],Ei__,2);
step_ = -p(2)/(2*p(1));
end