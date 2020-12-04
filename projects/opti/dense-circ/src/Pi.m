function Pi_ = Pi(r,lam,expo,i_)

nr = size(r,1);
jr_= 1:nr;
jr_(i_) = [];
Pi_ = 0;

for jr=1:(nr-1)
 j_ = jr_(jr);
 Pi__= (lam/sum( (r(i_,:)-r(j_,:)).^2 ))^expo;
 Pi_ = Pi_ + Pi__;
end

end