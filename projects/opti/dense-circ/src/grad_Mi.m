function g = grad_Mi(r,m,i_)

% nr = size(r,1);
% jr_= 1:nr;
% jr_(i_) = [];
% g = zeros(1,2);
% 
% for jr=1:(nr-1)
%  j_ = jr_(jr);
%  dist_ = sqrt( sum( (r(i_,:)-r(j_,:)).^2 ) );
%  li = 2*( dist_-d )/dist_;
%  g = g + li*(r(i_,:)-r(j_,:));
% end

dist_ = sqrt( sum( (r(i_,:)-m).^2 ) );
li = 1/dist_;
g = li*(r(i_,:)-m);
end