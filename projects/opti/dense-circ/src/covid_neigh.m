function [r,R,E,steps_] = covid_neigh(r,d,niter)

nr= size(r,1);
iter_= 0;

m = sum(r,1)/nr;
R = zeros(nr,2,niter);
E = zeros(niter,nr);
steps_=zeros(nr,2);

m = geome_median(r);

while iter_ < niter
 r_=r;
 
 % ------------------------------
 % minimize ||ri-m|| 
 % ------------------------------
 updates = zeros(nr,2);
 for ir=1:nr
  g = grad_Mi(r,m,ir);
  g = 0.25*g; % 0.5
 
  Mi_ = Mi(r,d,ir);
  k=[1e-4; 1e-3; 1e-2; 1e-1; 5e-1; 7.5e-1];
  step_ = step_Mi(r,g,m,Mi_,k,ir);
  if step_<0
    step_=0;
  end
  updates(ir,:) = -step_*g;
 end
 r=r+updates;
 
 % ------------------------------
 % minimize (lam/||ri-rj||^2)^m
 % ------------------------------
 lam=1;
 expo=2;
 for ir=1:nr
  g = grad_Pi(r,lam,expo,ir);
  g = 0.5*g; % 1
  
  Pi_ = Pi(r,lam,expo,ir);
  k=[1e-4; 1e-3; 1e-2; 1e-1];
  step_ = step_Pi(r,g,lam,expo,Pi_,k,ir);
  if step_<0
   step_=0;
  end
  E(iter_+1,ir) = Pi_;
  steps_(iter_+1,ir) = step_;
  r(ir,:) = r(ir,:) - step_*g;
  % updates(ir,:) = -step_*g;
 end
 
 % ------------------------------
 % minimize (||ri-rj||-d)^2
 % ------------------------------
 for ir=1:nr
  g = grad_Ei(r,d,ir);
  Ei_ = Ei(r,d,ir);
  k=[1e-4; 1e-3; 1e-2; 1e-1];
  step_ = step_Ei(r,g,d,Ei_,k,ir);
  if step_<0
    step_=0;
  end
  E(iter_+1,ir) = Ei_;
  steps_(iter_+1,ir) = step_;
  
  r(ir,:) = r(ir,:) - step_*g;
 end
 
 if or(sum(isnan(r(:)))>0 , sum(isinf(r(:)))>0)
  r=r_;
  return;
 end
 
 iter_=iter_+1;
 R(:,:,iter_) = r;
 
 if mod(iter_,fix(niter*0.05))==0
  fprintf('\n just finished iteration # %i',iter_);
 end
end
fprintf('\n\nbye bye\n\n')
end