function E = marko_E(w,r,rho_,miu_,betas,S)

ns= numel(w);
e_= ones(ns,1);

E = w.'*S*w + betas(1)*(r.'*w-rho_)^2 + betas(2)*(e_.'*w-miu_)^2;
end