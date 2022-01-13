function step_ = step_dista_(np,nedges,xy,edges,distao,gra_,iter,OBJ)
% diego domenzain 2022
nsteps= 1e2;
steps = logspace(-10,-1,nsteps);
objfnc= zeros(nsteps,1);

for isteps=1:nsteps
  % 👉
  dista = fwd_dista_(np,nedges,xy - steps(isteps)*gra_,edges);
  % 👇
  [objfnc_,err_]= obj_dista(dista, distao, OBJ);
  objfnc(isteps) = objfnc_;
end

[~,istep] = min(objfnc);
step_ = steps(istep);

% % 🐛
% figure(10);
% plot(steps,objfnc,'.-','markersize',20)
% title(num2str(iter))

end
