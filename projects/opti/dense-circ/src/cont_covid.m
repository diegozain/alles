% ------------------------------------------------------------------------------
r = covid_neigh(r,d);
% ------------------------------------------------------------------------------
D = dist_mat(r);
fprintf('\n max and min distances: %2.2d , %2.2d\n\n',max(D(:)),min( D(D(:)>0) ))
% ------------------------------------------------------------------------------
figure;
hold on
plot(r(:,1),r(:,2),'r.','markersize',50)
plot(r_(:,1),r_(:,2),'k.','markersize',20)
hold off
axis tight
axis image
xlabel('Length (m)')
ylabel('Width (m)')
title('Original vs Optimal')
simple_figure()
% ------------------------------------------------------------------------------
