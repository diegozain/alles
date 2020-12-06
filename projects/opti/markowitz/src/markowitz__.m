function [w,betas] = markowitz__(d,rho_,miu_)
S = cov(d);
r = mean(d,1);
r = r.';
ns= size(d,2);
e_= ones(ns,1);
A = [2*S -r -e_; r.' 0 0; e_.' 0 0];
b = [zeros(ns,1); rho_; miu_];
b = A\b;
w = b(1:ns);
betas=b((ns+1):(ns+2));
end