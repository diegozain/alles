function [w,W,E,steps_] = markowitz_(d,rho_,miu_,niter)

S = cov(d);
r = mean(d,1);
r = r.';
ns= size(d,2);
e_= ones(ns,1);
iter_=0;

W = [];
E = [];
steps_=[];

w = zeros(ns,1);

while iter_ < niter
 w_=w;
 % --------------------
 % sum to one 
 % --------------------
 g = 2*(e_.'*w -miu_)*e_;
 E_ = E_sum2one(w,miu_);
 k=[1e-4; 1e-3; 1e-2; 1e-1];
 step_ = step_sum2one(g,k,w,miu_,E_);
 if step_<0
   step_=0;
 elseif isnan(step_) 
   step_=1e-8;
 end
 dw= -step_*g;
 w = w + dw;
 % --------------------
 % return 
 % --------------------
 g = 2*(r.'*w - rho_)*r;
 E_ = E_return(w,r,rho_);
 k=[1e-4; 1e-3; 1e-2; 1e-1];
 step_ = step_return(g,k,w,r,rho_,E_);
 if step_<0
   step_=0;
 elseif isnan(step_) 
   step_=1e-8;
 end
 dw= -step_*g;
 w = w + dw;
 % --------------------
 % risk 
 % --------------------
 g = 2*S*w;
 E_ = E_risk(w,S);
 k=[1e-4; 1e-3; 1e-2; 1e-1];
 step_ = step_risk(g,k,w,S,E_);
 if step_<0
   step_=0;
 elseif isnan(step_) 
   step_=1e-8;
 end
 dw= -step_*g;
 w = w + dw;
 
 if or(sum(isnan(w))>0 , sum(isinf(w))>0)
  w=w_;
  fprintf('\n  your inversion went NaN (or Inf) at iteration # %i\n\n',iter_);
  return;
 end
 
 E_ = marko_E(w,r,rho_,miu_,[1,1],S);
 
 W(:,iter_+1) = w;
 E(iter_+1) = E_;
 steps_(iter_+1) = step_;
 
 W = [W w];
 E = [E E_];
 steps_ = [steps_ step_];
 
 iter_=iter_+1;
end

end