function step_ab = hd_step_ab(uo,uh,g_alphas,g_betas,Ob,k__,k___,nparabo,type_obj,t,alphas,betas,fos,h,nt_,nt__)
% ------------------------------------------------------------------------------
% --- parabola
% build perturbations
k_ = linspace(k__,k___,nparabo).';
Ob_ = zeros(nparabo+1,1);
Ob_(1) = Ob;
% compute many objective function values
for iparabo=1:nparabo
  % perturb
  alphas_= alphas - k_(iparabo)*g_alphas;
  betas_ = betas  - k_(iparabo)*g_betas;
  % fwd
  uh_ = hd_fwd(t,alphas_,betas_,fos,h,nt_,nt__);
  % obj
  [Ob,error_] = hd_obj(uo,uh_,type_obj);

  Ob_(iparabo+1) = Ob;
end
k_all = [0;k_];
% parabola approx
warning off;
p_ = polyfit(k_all,Ob_,2);
warning on;
% find zero of parabola (update = -step*gradient)
step_ab = -p_(2)/(2*p_(1));
end
