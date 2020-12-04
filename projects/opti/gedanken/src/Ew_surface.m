function EE = Ew_surface(p1,p2,s,M,d_o,N)
% ..............................................................................
np1 = length(p1);
np2 = length(p2);
EE = zeros(np1,np2);
% ..............................................................................
for i=1:np1
  for j=1:np2
    % point in parameter space
    %
    p = [p1(i) ; p2(j)];
    % forward model
    %
    [u,e] = fwd_w(p,s,M,d_o);
    % store in matrix
    %
    EE(i,j) = Ew_obj(e,N);
  end
end
% ..............................................................................
end
