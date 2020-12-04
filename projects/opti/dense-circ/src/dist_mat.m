function D = dist_mat(r)
nr = size(r,1);
D = zeros(nr,nr);
for i_=1:nr
 jr_= 1:nr;
 jr_(i_) = [];
 for jr=1:(nr-1)
  j_ = jr_(jr);
  D(i_,j_) = sqrt( sum( (r(i_,:)-r(j_,:)).^2 ) );
 end
end
end