function H = Hess_Ei(r,d,i_)

nr = size(r,1);
jr_= 1:nr;
jr_(i_) = [];
H = zeros(2,2);

for jr=1:(nr-1)
 j_ = jr_(jr);
 
 dist_ = sqrt( sum( (r(i_,:)-r(j_,:)).^2 ) );
 li = 2*( dist_-d )/dist_;
  
 dxli = (d*(r(i_,1) - r(j_,1)))/((r(i_,1) - r(j_,1))^2 + (r(i_,2) - r(j_,2))^2)^(3/2);
 dxli = 2*dxli;
 
 dyli = (d*(r(i_,2) - r(j_,2)))/((r(i_,1) - r(j_,1))^2 + (r(i_,2) - r(j_,2))^2)^(3/2);
 dyli = 2*dyli;
 
 H = H + [dxli*(r(i_,:)-r(j_,:)) ; dyli*(r(i_,:)-r(j_,:))] + [li 0 ; 0 li];
end
H=H.';
end
% ------------------------------------------------------------------------------
% compute dxli & dyli in Wolframalpha.com
% 
% [sqrt[[rix-rjx]^2 + [riy-rjy]^2] - d]/[ sqrt[[rix-rjx]^2 + [riy-rjy]^2] ]
% [sqrt[[a-b]^2     + [c-d]^2]     - d]/[ sqrt[[a-b]^2     + [c-d]^2] ]
% 
% dxli:
% 
% Diff[[sqrt[[a-b]^2 + [c-d]^2] - d]/[ sqrt[[a-b]^2 + [c-d]^2] ],a]
% (d (a - b))/((a - b)^2 + (c - d)^2)^(3/2)
% (d (rix - rjx))/((rix - rjx)^2 + (riy - rjy)^2)^(3/2)
% (d (r(i_,1) - r(j_,1)))/((r(i_,1) - r(j_,1))^2 + (r(i_,2) - r(j_,2))^2)^(3/2)
% 
% dyli:
% Diff[[sqrt[[a-b]^2 + [c-d]^2] - d]/[ sqrt[[a-b]^2 + [c-d]^2] ],c]
% (d (c - d))/((a - b)^2 + (c - d)^2)^(3/2)
% (d (riy - rjy))/((rix - rjx)^2 + (riy - rjy)^2)^(3/2)
% (d (r(i_,2) - r(j_,2)))/((r(i_,1) - r(j_,1))^2 + (r(i_,2) - r(j_,2))^2)^(3/2)
% ------------------------------------------------------------------------------


