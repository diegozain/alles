function Ei_ = Ei(r,d,i_)

nr = size(r,1);
jr_= 1:nr;
jr_(i_) = [];
Ei_ = 0;

for jr=1:(nr-1)
 j_ = jr_(jr);
 Ei__= (sqrt( sum( (r(i_,:)-r(j_,:)).^2 ) ) - d)^2;
 Ei_ = Ei_ + Ei__;
end

end