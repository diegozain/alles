function step_ = stepemcali(obj,param,datao,grad_,s,kparam_,kparam__,nparabo)
% diego domenzain
% jun 2022
% ------------------------------------------------------------------------------
% % --- pica. pica doesnt work for this inverse problem :(
% param_ = param + 1e-1*grad_;
% data_ = fwdemcali(param_,s);
% step_ = 1e-1 * ((data_.' * resi) / (data_.' * data_));
% ------------------------------------------------------------------------------
% --- build perturbations
% kparam = linspace(kparam_,kparam__,nparabo).';
kparam = logspace(log10(kparam_),log10(kparam__),nparabo).';
obj_ = zeros(nparabo+1,1);
obj_(1) = obj;
% compute many objective function values
for iparabo=1:nparabo
  % 🍐
  param_ = param - kparam(iparabo)*grad_;
  % 👉
  data_ = fwdemcali(param_,s);
  % Θ
  [obj,resi] = objemcali(data_,datao);
  obj_(iparabo+1) = obj;
end
kparam = [0;kparam];

if nparabo < 10
  % -- parabola approx
  warning off;
  p_ = polyfit(kparam,obj_,2);
  warning on;
  % find zero of parabola (update = -step*gradient)
  step_ = -p_(2)/(2*p_(1));
else
  % -- brute line-search
  [~,istep] = min(obj_);
  step_ = kparam(istep);
end
end
