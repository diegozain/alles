function E = E_sum2one(w,miu_)

ns= numel(w);
e_= ones(ns,1);

E = (e_.'*w-miu_)^2;
end