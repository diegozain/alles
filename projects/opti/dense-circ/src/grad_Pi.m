function g = grad_Pi(r,lam,expo,i_)

nr = size(r,1);
jr_= 1:nr;
jr_(i_) = [];
g = zeros(1,2);

for jr=1:(nr-1)
 j_ = jr_(jr);
 dist_ = sum( (r(i_,:)-r(j_,:)).^2 );
 li = -2*expo*lam*(lam/dist_)^(expo-1);
 g = g + li*(r(i_,:)-r(j_,:));
end
end