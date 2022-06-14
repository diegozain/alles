function [obj,resi] = objemcali(data,datao)
% diego domenzain
% jun 2022
% ------------------------------------------------------------------------------

% resi = data-datao;
% obj  = resi'*resi;

resi= data-datao;
obj = log(resi'*resi);
resi= (1/obj) * resi;

end
