function [alphas,betas,fos] = hd_inversion(uo,niter,t,alphas,betas,fos,h,nt_,nt__)
% ------------------------------------------------------------------------------
% diego domenzain. jul 2021, @ AU.
% ------------------------------------------------------------------------------
% uo     : observed data with harmonic noise
% niter  : # of iterations
% t      : time vector
% alphas : alpha coefficients for harmonic series. size (nb*nh x 1)
% betas  : beta coefficients for harmonic series.  size (nb*nh x 1)
% fos    : frequencies that are harmonized.        size (nb x 1)
% h      : multiples of fo that make up harmonics. size (1 x nh)
% nt_    : # of samples in the over-lapping time
% nt__   : # of overlapping time samples   (nt__ = ceil(1/fo/dt))
% ------------------------------------------------------------------------------
% nb     : # of time blocks
% nh     : # of harmonics
% ------------------------------------------------------------------------------
nb = numel(fos);
nh = numel(h);
% ------------------------------------------------------------------------------
% -- parabola step sizes hyperparameters

% - synthetic example
% fo
k_fos_     =1e-8;
k_fos__    =2e-5; % 1e-6
nparabo_fos=1e2;
% α & β
k_alphas_ =1e-8;
k_alphas__=1e-3;
k_betas_  =1e-8;
k_betas__ =5e-3;
nparabo_ab=1e2;

% % - field data
% % fo
% k_fos_    =1e-9;
% k_fos__   =1e-8; % 1e-6
% nparabo_fos=1e2;
% % α & β
% k_alphas_ =1e-8;
% k_alphas__=1e-2;% 1e-3;
% k_betas_  =1e-8;
% k_betas__ =1e-2;% 5e-3;
% nparabo_ab=1e2;
% ------------------------------------------------------------------------------
% -- memory over iterations
% fo
ob_fos   =zeros(niter,1);
fos_niter=zeros(nb,niter);
% α & β
ob_alphas=zeros(niter,1);
alphas_niter=zeros(nb*nh,niter);
ob_betas=zeros(niter,1);
betas_niter=zeros(nb*nh,niter);
% ------------------------------------------------------------------------------
for iiter=1:niter
  % -- fos
  % type of objective function.
  % 'sse'   = sum of squared errors
  % 'lnsse' = ln( sum of squared errors )
  type_obj='lnsse';
  uh          = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
  [Ob,error_] = hd_obj(uo,uh,type_obj);

  ob_fos(iiter,1)=Ob;
  fos_niter(:,iiter)=fos;
  % fos gradient
  g_fos = hd_grad_f(error_,t,alphas,betas,fos,h,nt_,nt__);
  % fos step size
  step_fos = hd_step_f(uo,uh,g_fos,Ob,k_fos_,k_fos__,nparabo_fos,type_obj,t,alphas,betas,fos,h,nt_,nt__);
  % update
  d_fos = - step_fos * g_fos;
  % frequency always positive
  fos   = fos.*exp(fos.*d_fos);
end
[~,iiter] = min(ob_fos);
fos = fos_niter(:,iiter);

for iiter=1:niter
  % type of objective function.
  % 'sse'   = sum of squared errors
  % 'lnsse' = ln( sum of squared errors )
  type_obj='sse';
  uh          = hd_fwd(t,alphas,betas,fos,h,nt_,nt__);
  [Ob,error_] = hd_obj(uo,uh,type_obj);

  ob_alphas(iiter,1)=Ob;
  alphas_niter(:,iiter)=alphas;

  ob_betas(iiter,1)=Ob;
  betas_niter(:,iiter)=betas;

  % -- alphas
  % alphas gradient
  g_alphas = hd_grad_a(error_,t,alphas,fos,h,nt_,nt__);
  % alphas step size
  step_alphas = hd_step_a(uo,uh,g_alphas,Ob,k_alphas_,k_alphas__,nparabo_ab,type_obj,t,alphas,betas,fos,h,nt_,nt__);

  % -- betas
  % betas gradient
  g_betas = hd_grad_b(error_,t,betas,fos,h,nt_,nt__);
  % betas step size
  step_betas = hd_step_b(uo,uh,g_betas,Ob,k_betas_,k_betas__,nparabo_ab,type_obj,t,alphas,betas,fos,h,nt_,nt__);

  % -- update
  % alphas
  d_alphas = - step_alphas * g_alphas;
  alphas   = alphas + d_alphas;
  % betas
  d_betas  = - step_betas * g_betas;
  betas    = betas + d_betas;
end
[~,iiter] = min(ob_alphas);
alphas = alphas_niter(:,iiter);
[~,iiter] = min(ob_betas);
betas = betas_niter(:,iiter);
end