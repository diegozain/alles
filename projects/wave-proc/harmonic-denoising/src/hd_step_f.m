function step_fos = hd_step_f(uo,uh,g_fos,Ob,k_fos_,k_fos__,nparabo,type_obj,t,alphas,betas,fos,h,nt_,nt__)
% ------------------------------------------------------------------------------
% --- build perturbations
% k_fos = linspace(k_fos_,k_fos__,nparabo).';
k_fos = logspace(log10(k_fos_),log10(k_fos__),nparabo).';
Ob_ = zeros(nparabo+1,1);
Ob_(1) = Ob;
% compute many objective function values
for iparabo=1:nparabo
 % perturb
 fos_ = fos.*exp(- k_fos(iparabo)*g_fos.*fos);
 % fwd
 uh_  = hd_fwd(t,alphas,betas,fos_,h,nt_,nt__);
 % obj
 [Ob,error_] = hd_obj(uo,uh_,type_obj);

 Ob_(iparabo+1) = Ob;
end
k_fos = [0;k_fos];

if nparabo < 10
  % -- parabola approx
  warning off;
  p_ = polyfit(k_fos,Ob_,2);
  warning on;
  % find zero of parabola (update = -step*gradient)
  step_fos = -p_(2)/(2*p_(1));
else
  % -- brute line-search
  [~,istep] = min(Ob_);
  step_fos = k_fos(istep);
end
% % -- for debugging
% figure;
% hold on;
% plot(k_fos,Ob_,'r.-','markersize',30);
% % plot(linspace(min(k_fos),max(k_fos),1e3),polyval(p_,linspace(min(k_fos),max(k_fos),1e3)))
% % hold off;
end
