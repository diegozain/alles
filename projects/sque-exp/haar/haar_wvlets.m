function [phi_ij,phi_ij_] = haar_wvlets(t,ii,jj)
  % ii *has* to be any number of the form:
  % ii = 0, 1, ..., (2^jj - 1)
  % 
  % jj can be any positive integer.
 
 nt = numel(t);
 
 phi_ij = haar_phi((2^jj)*t-ii,nt);
 phi_ij = (2^(jj*0.5))*phi_ij;
 
 phi_ij_= haar_phi_((2^jj)*t-ii,nt);
 phi_ij_= (2^(jj*0.5))*phi_ij_;
 
end
% ------------------------------------------------------------------------------

function phi = haar_phi(t,nt)
  phi= zeros(nt,1);
  for it=1:nt
    a=0;
    if 0 < t(it) && t(it) <= 1
      a=1;
    end
    phi(it) = a;
  end
end

% ------------------------------------------------------------------------------

function phi_ = haar_phi_(t,nt)
 phi_= zeros(nt,1);
 for it=1:nt
   a = 0;
   if 0 < t(it) && t(it) <= 0.5
     a=1;
   elseif 0.5 < t(it) && t(it) <= 1
     a=-1;
   end
   phi_(it) = a;
 end
end

% ------------------------------------------------------------------------------
