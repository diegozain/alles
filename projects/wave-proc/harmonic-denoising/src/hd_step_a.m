function step_alphas = hd_step_a(uo,uh,g_alphas,Ob,k_alphas_,k_alphas__,nparabo,type_obj,t,alphas,betas,fos,h,nt_,nt__)
% % --- pica. pica doesnt work for this inverse problem :( 
% alphas_     = alphas + k_alphas*g_alphas;
% uh_         = hd_fwd(t,alphas_,betas,fos,h,nt_,nt__);
% step_alphas = k_alphas * ((uh_.' * error_) / (uh_.' * uh_));
% ------------------------------------------------------------------------------
% --- build perturbations
k_alphas = linspace(k_alphas_,k_alphas__,nparabo).';
Ob_ = zeros(nparabo+1,1);
Ob_(1) = Ob;
% compute many objective function values
for iparabo=1:nparabo
  % perturb
  alphas_ = alphas - k_alphas(iparabo)*g_alphas;
  % fwd
  uh_ = hd_fwd(t,alphas_,betas,fos,h,nt_,nt__);
  % obj
  [Ob,error_] = hd_obj(uo,uh_,type_obj);
  
  Ob_(iparabo+1) = Ob;
end
k_alphas = [0;k_alphas];

if nparabo < 10
  % -- parabola approx
  warning off;
  p_ = polyfit(k_alphas,Ob_,2);
  warning on;
  % find zero of parabola (update = -step*gradient)
  step_alphas = -p_(2)/(2*p_(1));
else
  % -- brute line-search
  [~,istep] = min(Ob_);
  step_alphas = k_alphas(istep);
end
% figure;
% plot(k_alphas,Ob_,'r.-','markersize',30);
end