function step_ = step_dista(x,r,g,iter,OBJ)
% diego domenzain 2021
nsteps= 1e2;
steps = logspace(-10,-1,nsteps);
objfnc= zeros(nsteps,1);

for isteps=1:nsteps
  % 👉
  dat_ = fwd_dista(x - steps(isteps)*g);
  % 👇
  [objfnc_,err_]= obj_dista(dat_, r, OBJ);
  objfnc(isteps) = objfnc_;
end

[~,istep] = min(objfnc);
step_ = steps(istep);

% % 🐛
% figure(10);
% plot(steps,objfnc,'.-','markersize',20)
% title(num2str(iter))

end
