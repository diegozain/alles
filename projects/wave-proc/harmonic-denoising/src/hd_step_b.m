function step_betas = hd_step_b(uo,uh,g_betas,Ob,k_betas_,k_betas__,nparabo,type_obj,t,alphas,betas,fos,h,nt_,nt__)
% % --- pica. pica doesnt work for this inverse problem :(
% % g_betas    = normali(g_betas);
% betas_     = betas + k_betas*g_betas;
% uh_        = hd_fwd(t,alphas,betas_,fos,h,nt_,nt__);
% step_betas = k_betas * ((uh_.' * error_) / (uh_.' * uh_));
% ------------------------------------------------------------------------------
% --- parabola
% build perturbations
k_betas = linspace(k_betas_,k_betas__,nparabo).';
Ob_ = zeros(nparabo+1,1);
Ob_(1) = Ob;
% compute many objective function values
for iparabo=1:nparabo
  % perturb
  betas_ = betas - k_betas(iparabo)*g_betas;
  % fwd
  uh_ = hd_fwd(t,betas_,betas,fos,h,nt_,nt__);
  % obj
  [Ob,error_] = hd_obj(uo,uh_,type_obj);
  
  Ob_(iparabo+1) = Ob;
end
k_betas = [0;k_betas];

if nparabo < 10
  % -- parabola approx
  warning off;
  p_ = polyfit(k_betas,Ob_,2);
  warning on;
  % find zero of parabola (update = -step*gradient)
  step_betas = -p_(2)/(2*p_(1));
else
  % -- brute line-search
  [~,istep] = min(Ob_);
  step_betas = k_betas(istep);
end
% figure;
% plot(sort(k_betas),Ob_,'r.-','markersize',30);
end