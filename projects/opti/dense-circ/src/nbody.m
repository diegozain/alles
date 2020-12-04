function [P,M] = nbody(P,M,dt,nt)

np = size(P,1);
r  = P(:,:,1);
v  = P(:,:,2);
a  = zeros(np,2);

nei= 1:np;
% figure;
for it=1:nt
  for ip=1:np
    nei_ = nei;
    nei_(ip) = [];
    for jp=1:(np-1)
     jp_ = nei_(jp);
     dist_ = sum((r(ip,:) - r(jp_,:)).^2);
     a(ip,:) = a(ip,:) - ((M(jp_)/ dist_ ) * (r(ip,:) - r(jp_,:)));
    end
  end
  
  for ip=1:np
    r(ip,:) = r(ip,:) + v(ip,:)*dt;
    v(ip,:) =  a(ip,:)*dt;
  end
  
  % plot(r(:,1),r(:,2),'k.','markersize',50)
  % xlim([-5 5])
  % ylim([-5 8])
  % pause
 
  a=zeros(np,2);
end

P(:,:,1) = r;
P(:,:,2) = v;

end